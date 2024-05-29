﻿//------------------------------------------------------------------------
// Copyright(c) 2024 yg331.
//------------------------------------------------------------------------

#pragma once
#include "comp_cids.h"
#include "comp_dataexchange.h"

#include "public.sdk/source/vst/vstaudioeffect.h"
#include "pluginterfaces/base/ustring.h"
#include "base/source/fstring.h"

#include <algorithm>
#include <limits>
#include <cmath>
#include <vector>
#include <numeric>
#include <queue>
#include <deque>
#ifndef M_PI
#define M_PI        3.14159265358979323846264338327950288   /* pi             */
#endif

namespace yg331 {


    class Kaiser {
    public:

#define maxTap 512

        static inline double Ino(double x) 
        {
            double d = 0, ds = 1, s = 1;
            do
            {
                d += 2;
                ds *= x * x / (d * d);
                s += ds;
            } while (ds > s * 1e-6);
            return s;
        };

        static void calcFilter(double Fs, double Fa, double Fb, int M, double Att, double* dest) 
        {
            // Kaiser windowed FIR filter "DIGITAL SIGNAL PROCESSING, II" IEEE Press pp 123-126.

            int Np = (M - 1) / 2;
            double A[maxTap] = { 0, };
            double Alpha;
            double Inoalpha;

            A[0] = 2 * (Fb - Fa) / Fs;

            for (int j = 1; j <= Np; j++)
                A[j] = (sin(2.0 * j * M_PI * Fb / Fs) - sin(2.0 * j * M_PI * Fa / Fs)) / (j * M_PI);

            if (Att < 21.0)
                Alpha = 0;
            else if (Att > 50.0)
                Alpha = 0.1102 * (Att - 8.7);
            else
                Alpha = 0.5842 * pow((Att - 21), 0.4) + 0.07886 * (Att - 21);

            Inoalpha = Ino(Alpha);

            for (int j = 0; j <= Np; j++)
            {
                dest[Np + j] = A[j] * Ino(Alpha * std::sqrt(1.0 - ((double)(j * j) / (double)(Np * Np)))) / Inoalpha;
            }
            for (int j = 0; j < Np; j++)
            {
                dest[j] = dest[M - 1 - j];
            }

        };
    };

    class Flt {
    public:
        double coef alignas(16)[maxTap] = { 0, };
        double buff alignas(16)[maxTap] = { 0, };
        int now = 0;
        int size = maxTap;
        double* buff_ptr = buff;
        void acc() {
            now++;
            if (now == size) {
                now = 0;
            }
        }
        int get_nth(int n) {
            if (now + n >= size) {
                return now + n - size;
            }
            else {
                return now + n;
            }
        }
    };

    class smooth {
    public: 
        static double ABS(double in, double smoothness) 
        {
            double err = sqrt(1.0 + smoothness * smoothness) - 1.0;
            return sqrt((in * in) + smoothness * smoothness) - err; //sqrt(smoothness);
        };

        static double Max(double a, double b, double smoothness) 
        {
            return (a + b + ABS(a - b, smoothness)) * 0.5;
        };

        static double Step(double in, double smoothness, double epsilon = 0.0)
        {
            double init = 1.0 / (1.0 + exp(-1.0 * smoothness * in));
            if (1.0 - init < epsilon) init = 1.0;
            else if (init < epsilon) init = 0.0;
            return init;
        };
    };

    class Decibels
    {
    public:
        //==============================================================================
        /** Converts a dBFS value to its equivalent gain level.

            A gain of 1.0 = 0 dB, and lower gains map onto negative decibel values. Any
            decibel value lower than minusInfinityDb will return a gain of 0.
        */
        template <typename Type>
        static Type decibelsToGain(Type decibels,
            Type minusInfinityDb = Type(defaultMinusInfinitydB))
        {
            return decibels > minusInfinityDb ? std::pow(Type(10.0), decibels * Type(0.05))
                : Type();
        }

        /** Converts a gain level into a dBFS value.

            A gain of 1.0 = 0 dB, and lower gains map onto negative decibel values.
            If the gain is 0 (or negative), then the method will return the value
            provided as minusInfinityDb.
        */
        template <typename Type>
        static Type gainToDecibels(Type gain,
            Type minusInfinityDb = Type(defaultMinusInfinitydB))
        {
            return gain > Type() ? std::max(minusInfinityDb, static_cast<Type> (std::log10(gain)) * Type(20.0))
                : minusInfinityDb;
        }

        //==============================================================================
        /** Converts a decibel reading to a string.

            By default the returned string will have the 'dB' suffix added, but this can be removed by
            setting the shouldIncludeSuffix argument to false. If a customMinusInfinityString argument
            is provided this will be returned if the value is lower than minusInfinityDb, otherwise
            the return value will be "-INF".
        */
        

    private:
        //==============================================================================
        enum { defaultMinusInfinitydB = -100 };

        Decibels() = delete; // This class can't be instantiated, it's just a holder for static methods..
    };



    /*CrestFactor Class:
    * Calculates the average Crest-Factor for a given buffer.
    * Crest-Factor is time-variable value calculated from the ratio between peak and rms of the signal
    */
    class CrestFactor
    {
    public:

        CrestFactor() = default;

        // Prepares processor with ProcessSpec-Object and recalculates coefficients for current ballistics
        void prepare(const double& fs) {
            sampleRate = fs;
            //Calculate alpha for release time of 200ms, same release time for peak & rms detector
            a1 = exp(-1.0 / (sampleRate * 0.2));
            b1 = 1 - a1;
        };

        // Calculates Crest-Factor for given buffer
        void process(const double* src, int numSamples) {
            //Init accumulators
            if (!peakState) peakState = src[0];
            if (!rmsState) rmsState = src[0];

            //Reset avg attack/release
            avgAttackTime = 0.0;
            avgReleaseTime = 0.0;

            //Calculate averages of auto - attack/release times for a single buffer
            for (int i = 0; i < numSamples; i++)
            {
                //Square of input signal
                const double s = static_cast<double>(src[i]) * static_cast<double>(src[i]);

                //Update peak state
                peakState = (std::max)(s, a1 * peakState + b1 * s);

                //Update rms state
                rmsState = a1 * rmsState + b1 * s;

                //calculate squared crest factor
                const double c = peakState / rmsState;
                cFactor = c > 0.0 ? c : 0.0;

                //calculate ballistics
                if (cFactor > 0.0)
                {
                    attackTimeInSeconds = 2 * (maxAttackTime / cFactor);
                    releaseTimeInSeconds = 2 * (maxReleaseTime / cFactor) - attackTimeInSeconds;

                    //Update avg ballistics
                    avgAttackTime += attackTimeInSeconds;
                    avgReleaseTime += releaseTimeInSeconds;
                }
            }

            // Calculate average ballistics & crest factor
            avgAttackTime /= numSamples;
            avgReleaseTime /= numSamples;
        };

        // Get average calculated attack time of a buffer, call after proces()
        double getAvgAttack() { return avgAttackTime; };

        // Get average calculated release time of a buffer, call after process()
        double getAvgRelease() { return avgReleaseTime; };

    private:
        double attackTimeInSeconds{ 0.0 };
        double releaseTimeInSeconds{ 0.14 };
        double avgAttackTime{ 0.0 };
        double avgReleaseTime{ 0.14 };

        double peakState{ 0.0 };
        double rmsState{ 0.0 };
        double a1{ 0.0 }, b1{ 0.0 };
        double sampleRate{ 0.0 };
        double maxAttackTime{ 0.08 }, maxReleaseTime{ 1.0 }; //respective 8ms and 1sec
        double cFactor{ 0.0 };
    };


    /*Simple exponential moving average filter, also known as 1-pole iir filter
    * This class can be used to smooth values over a certain time frame
    */
    class SmoothingFilter
    {
    public:

        SmoothingFilter() = default;

        // Prepares the SmoothingFilter with a sampleRate
        void prepare(const double& fs) {
            sampleRate = fs;
            a1 = 1;
            b1 = 1 - a1;
        };

        // Processes a given sample
        void process(const double& sample) {
            if (first)
            {
                state = sample;
                first = false;
            }
            state = a1 * sample + b1 * state;
        };

        // Sets coefficient manually
        void setAlpha(double a) {
            a1 = a;
            b1 = 1 - a1;
        };

        // Set time-frame in seconds, recalculates needed coefficients
        void setAlphaWithTime(double timeInSeconds) {
            a1 = exp(-1.0 / (sampleRate * timeInSeconds)); 
            b1 = 1 - a1;
        };

        // Gets current value
        double getState() {
            return state;
        };

    private:
        double a1{ 1.0 }, b1{ 0.0 };
        double state{ 0.0 };
        double sampleRate{ 0.0 };
        bool first{ true };
    };


    /*LevelDetector Class:
    * Used to have a smooth representation of the level Might be used in linear or log. 
    * Domain In this compressor implementation it's used in log. 
    * Domain after the gain computer to smooth the calculated attenuations,
    * Therefore the detector does not have to work on the whole dynamic range of the input signal
    */
    class LevelDetector
    {
    public:
        LevelDetector() = default;

        // Prepares LevelDetector with a ProcessSpec-Object containing samplerate, blocksize and number of channels
        void prepare(const double& fs) {
            sampleRate = fs;
            crestFactor.prepare(fs);
            attackSmoothingFilter.prepare(fs);
            releaseSmoothingFilter.prepare(fs);

            alphaAttack = exp(-1.0 / (sampleRate * attackTimeInSeconds));
            alphaRelease = exp(-1.0 / (sampleRate * releaseTimeInSeconds));
            state01 = 0.0;
            state02 = 0.0;
        };

        // Sets attack time constant
        void setAttack(const double& attack) {
            if (attack != attackTimeInSeconds)
            {
                attackTimeInSeconds = attack; //Time it takes to reach 1-1/e = 0.63
                alphaAttack = exp(-1.0 / (sampleRate * attackTimeInSeconds)); //aA = e^(-1/TA*fs)

                double w = 1 / (2.0 * attackTimeInSeconds * sampleRate);
                g_a = tan(w);
                gt0_a = 1 / (1 + g_a);
                gt1_a = g_a * gt0_a;
            }
        };

        // Sets release time constant
        void setRelease(const double& release) {
            if (release != releaseTimeInSeconds)
            {
                releaseTimeInSeconds = release; //Time it takes to reach 1 - (1-1/e) = 0.37
                alphaRelease = exp(-1.0 / (sampleRate * releaseTimeInSeconds)); //aR = e^(-1/TR*fs)

                // a1 = d, tau = (sampleRate * timeInSeconds), f = 1 / 2pi*tau
                //double tau = releaseTimeInSeconds;
                //double f = 1 / 2 * M_PI * tau;
                //double w = f * M_PI / sampleRate;
                // w = 2 * pi * f = 1 / tau
                double w = 1 / (2.0 * releaseTimeInSeconds * sampleRate);
                // g = tan (pi * f / sample rate)

                g_r = tan(w);
                gt0_r = 1 / (1 + g_r);
                gt1_r = g_r * gt0_r;
            }
        };

        // Sets auto attack to enabled/disabled
        void setAutoAttack(bool isEnabled) {  autoAttack = isEnabled; };

        // Sets auto release to enabled/disabled
        void setAutoRelease(bool isEnabled) { autoRelease = isEnabled; };

        // Gets current attack time constant
        double getAttack() {  return attackTimeInSeconds; };

        // Gets current release time constant
        double getRelease() {  return releaseTimeInSeconds; };

        // Gets calculated attack coefficient
        double getAlphaAttack() { return alphaAttack; };

        // gets calculated release coefficient
        double getAlphaRelease() { return alphaRelease; };

        // Processes a sample with smooth branched peak detector
        // 1. Update state with coef branch, instead of switching state.
        // 2. Smooth transition on switching branches. 
        double processPeakBranched(const double& in) {
            
            double d = state01 - in;
            double pp = smooth::Step(d, 100.0, 0.0); // (delta > 0) ? 1.0 : 0.0;
            pp *= pp;
            double nn = 1.0 - pp;
            state01 = (pp * alphaAttack + nn * alphaRelease) * state01 
                    + (1 - (pp * alphaAttack + nn * alphaRelease)) * in;
            return static_cast<double>(state01); //y_L
            

            /*
            double delta = v2 - in; 
            double p = smooth::Step(delta, 100.0); // (delta > 0) ? 1.0 : 0.0;
            p *= p;
            double n = 1.0 - p;
            v2 = (p * gt1_a + n * gt1_r) * in 
               + (p * gt0_a + n * gt0_r) * ic1eq;
            ic1eq += 2.0 * (p * g_a + n * g_r) * (in - v2);
            //v2 = gt1 * in + gt0 * ic1eq;
            //ic1eq += 2.0 * g * (in - v2);
            return static_cast<double>(v2); 
            */

            //Smooth branched peak detector
            if (in < state01)
                state01 = alphaAttack * state01 + (1 - alphaAttack) * in;
            else
                state01 = alphaRelease * state01 + (1 - alphaRelease) * in;

            return static_cast<double>(state01); //y_L
        };

        // Processes a sample with smooth decoupled peak detector
        double processPeakDecoupled(const double& in) {

            const double input = static_cast<double>(in);
            /*
            v1 = (std::max)(input, gt1_r * input + gt0_r * ic2eq);
            ic2eq += 2.0 * g_r * (input - v1);
            v2 = gt1_a * v1 + gt0_a * ic1eq;
            ic1eq += 2.0 * g_a * (v1 - v2);
            return static_cast<double>(v2); //y_L
            */
            //float smooth_max(float a, float b) { return (a + b + smooth_ABS(a - b)) * 0.5f; };
            //Smooth decoupled peak detector
            //const double input = static_cast<double>(in);
            state02 = (std::max)(input, alphaRelease * state02 + (1 - alphaRelease) * input);
            state01 = alphaAttack * state01 + (1 - alphaAttack) * state02;
            return static_cast<double>(state01);
        };

        // Applies ballistics to given buffer
        void applyBallistics(double* src, int numSamples) {
            // Apply ballistics to src buffer
            for (int i = 0; i < numSamples; i++) {
                //src[i] = processPeakDecoupled(src[i]);
                src[i] = processPeakBranched(src[i]);
            }
        };

        // Processes crest factor and sets ballistics accordingly
        void processCrestFactor(const double* src, int numSamples) {
            if (autoAttack || autoRelease)
            {
                //Crest factor calculation
                crestFactor.process(src, numSamples);
                attackSmoothingFilter.process(crestFactor.getAvgAttack());
                releaseSmoothingFilter.process(crestFactor.getAvgRelease());
                if (autoAttack) setAttack(attackSmoothingFilter.getState());
                if (autoRelease) setRelease(releaseSmoothingFilter.getState());
            }
        };

    private:
        CrestFactor crestFactor;
        SmoothingFilter attackSmoothingFilter;
        SmoothingFilter releaseSmoothingFilter;

        double attackTimeInSeconds{ 0.01 }, alphaAttack{ 0.0 };
        double releaseTimeInSeconds{ 0.14 }, alphaRelease{ 0.0 };
        double state01{ 0.0 }, state02{ 0.0 };
        double sampleRate{ 0.0 };
        bool autoAttack{ false };
        bool autoRelease{ false };

        double g_a = 0.5, g_r = 0.5;
        double gt0_a = 1 / (g_a), gt0_r = 1 / (g_r);
        double gt1_a = g_a * gt0_a;
        double gt1_r = g_r * gt0_r;
        double v1 = 0, ic1eq = 0;
        double v2 = 0, ic2eq = 0;
    };

    /* GainComputer Class:
     * Calculates the needed attenuation to compress a signal with given characteristics
     */
    class GainComputer
    {
    public:

        GainComputer() {
            threshold = -20.0;
            ratio = 2.0;
            slope = 1.0 / ratio - 1.0;
            knee = 6.0;
            kneeHalf = 3.0;
        };

        // Sets the threshold in dB
        void setThreshold(double newTreshold)
        {
            threshold = newTreshold;
        };

        // Sets the ratio in dB
        void setRatio(double newRatio)
        {
            if (ratio != newRatio)
            {
                ratio = newRatio;
                if (ratio > 23.9) ratio = -std::numeric_limits<double>::infinity();
                slope = 1.0 / newRatio - 1.0;
            }
        };

        // Sets the knee-width in dB (if > 0, 2nd order interpolation for soft knee)
        void setKnee(double newKnee)
        {
            if (newKnee != knee)
            {
                knee = newKnee;
                kneeHalf = newKnee / 2.0;
            }
        };

        // Applies characteristics to a given sample, 2nd order spline interpolation
        // returns attenuation Xl == Xg - Yg
        double applyCompression(double& input)
        {
            
            const double overshoot = input - threshold;

            if (overshoot <= -kneeHalf)
                return 0.0;
            if (overshoot > -kneeHalf && overshoot <= kneeHalf)
                return 0.5 * slope * ((overshoot + kneeHalf) * (overshoot + kneeHalf)) / knee;

            return slope * overshoot; 
        };

        void applyCompressionToBuffer(double* src, int numSamples)
        {
            for (int i = 0; i < numSamples; i++)
            {
                //const double level = std::max(abs(src[i]), 1e-6);
                double levelInDecibels = Decibels::gainToDecibels(src[i]);
                src[i] = applyCompression(levelInDecibels);
            }
        };

    private:
        double threshold{ -20.0f };
        double ratio{ 2.0f };
        double knee{ 6.0f }, kneeHalf{ 3.0f };
        double slope{ -0.5f };
    };

    
    /* Compressor-Class:
     * The circruit is modeled after the "ideal" VCA-Compressor
     * based on the paper "Digital Dynamic Range Compressor Design �  Tutorial and Analysis"
     * by Giannoulis, Massberg & Reiss
     */
    class Compressor
    {
    public:
        Compressor() = default;
        ~Compressor() {
            rawSidechainSignal = nullptr;
        };

        // Prepares compressor with a ProcessSpec-Object containing samplerate, blocksize and number of channels
        void prepare(const double& sampleRate, const double& maximumBlockSize) {
            _sampleRate = sampleRate;
            _maximumBlockSize = maximumBlockSize;
            ballistics.prepare(2 * sampleRate);
            originalSignal[0].resize(maximumBlockSize, 0.0);
            originalSignal[1].resize(maximumBlockSize, 0.0);
            sidechainSignal.resize(maximumBlockSize, 0.0);
            rawSidechainSignal = sidechainSignal.data();
            smoothedAutoMakeup.prepare(sampleRate);
            smoothedAutoMakeup.setAlpha(0.03);

            for (int ch = 0; ch < 2; ch++) {
                up[ch].size = fir_size;
                dn[ch].size = fir_size;
                Kaiser::calcFilter(96000.0, 0.0, 24000.0, fir_size, 180.0, up[ch].coef);
                Kaiser::calcFilter(96000.0, 0.0, 24000.0, fir_size, 180.0, dn[ch].coef);
                for (int i = 0; i < fir_size; i++) up[ch].coef[i] *= 2.0;
            }
            upsample_buff[0].resize(2 * maximumBlockSize, 0.0);
            upsample_buff[1].resize(2 * maximumBlockSize, 0.0);
            dnsample_buff[0].resize(2 * maximumBlockSize, 0.0);
            dnsample_buff[1].resize(2 * maximumBlockSize, 0.0);
        };

        // Sets compressor to bypassed/not bypassed
        void setPower(bool newPower) { bypassed = newPower; };

        // Sets input in dB
        void setInput(double newInput) { input = newInput; };

        // Sets threshold in dB
        void setThreshold(double thresholdInDb) { gainComputer.setThreshold(thresholdInDb); };

        // Sets ratio in dB
        void setRatio(double rat) { gainComputer.setRatio(rat); };

        // Sets knee-width in dB (> 0 = soft knee)
        void setKnee(double kneeInDb) { gainComputer.setKnee(kneeInDb); };


        void setLookAhead(double _LookAhead) { lookaheadDelay = _LookAhead; };
        void setLookAheadEnable(bool _LookAheadEnable) { LookAheadEnable = _LookAheadEnable; };

        // Sets make-up gain in dB
        void setMakeup(double makeupGainInDb) { makeup = makeupGainInDb; };

        // Sets mix 0.0f - 1.0f
        void setMix(double newMix) { mix = newMix; };

        // Sets attack time in milliseconds
        void setAttack(double attackTimeInMs) { ballistics.setAttack(attackTimeInMs * 0.001); };

        // Sets release time in milliseconds
        void setRelease(double releaseTimeInMs) { ballistics.setRelease(releaseTimeInMs * 0.001); };

        // Sets auto attack to enabled = true or disabled = false
        void setAutoAttack(bool isEnabled) { ballistics.setAutoAttack(isEnabled); };

        // Sets auto release to enabled = true or disabled = false
        void setAutoRelease(bool isEnabled) { ballistics.setAutoRelease(isEnabled); };

        // Sets auto makeup to enabled = true or disabled = false
        void setAutoMakeup(bool isEnabled) { autoMakeupEnabled = isEnabled; };

        // Gets current make-up gain value
        double getMakeup() { return makeup; };

        // Return current sampleRate
        double getSampleRate() { return _sampleRate; };

        double getMaxGainReduction() { return maxGainReduction; };

        // Processes input buffer
        template <typename SampleType>
        void process(
            SampleType** inputs,
            SampleType** outputs,
            Steinberg::int32 numChannels,
            Steinberg::Vst::SampleRate getSampleRate,
            Steinberg::int32 sampleFrames)
        {
            if (!bypassed)
            {
                // Clear any old samples
                for (int i = 0; i < sampleFrames; i++) {
                    originalSignal[0][0] = 0.0;
                    originalSignal[1][0] = 0.0;
                    rawSidechainSignal[i] = 0.0f;
                }
                maxGainReduction = 0.0f;

                // Apply input gain
                for (int ch = 0; ch < numChannels; ch++) applyInputGain(inputs[ch], sampleFrames);

                // Get max l/r amplitude values and fill sidechain signal
                for (int ch = 0; ch < numChannels; ch++) {
                    for (int i = 0; i < sampleFrames; i++) {
                        // Multi-phase all-pass interpolation
                        // https://www.dsprelated.com/freebooks/pasp/First_Order_Allpass_Interpolation.html
                        // p1_sy[ch] = 0.198912367379658 * inputs[ch][i] + p1_sx[ch] - 0.198912367379658 * p1_sy[ch]; // 15000/48000
                        // p1_sx[ch] = inputs[ch][i];
                        // p2_sy[ch] = 0.41421356237309503 * inputs[ch][i] + p2_sx[ch] - 0.41421356237309503 * p2_sy[ch]; // 18000/48000
                        // p2_sx[ch] = inputs[ch][i];
                        // p3_sy[ch] = 0.6681786379192988 * inputs[ch][i] + p3_sx[ch] - 0.6681786379192988 * p3_sy[ch]; // 21000/48000
                        // p3_sx[ch] = inputs[ch][i];

                        //double a = std::max(abs((double)inputs[ch][i]), abs(p1_sx[ch]));
                        //double b = std::max(abs(p2_sy[ch]), abs(p3_sy[ch]));
                        //double c = std::max(abs(p4_sy[ch]), abs(p5_sy[ch]));
                        //double d = std::max(a, std::max(b, c));
                        //double d = std::max(a, b);
                        
                        //rawSidechainSignal[i] = (rawSidechainSignal[i] > d) ? rawSidechainSignal[i] : d;
                        
                        
                        rawSidechainSignal[i] = (rawSidechainSignal[i] > abs((double)inputs[ch][i]))
                            ? rawSidechainSignal[i] : (double)inputs[ch][i];
                        //rawSidechainSignal[i] = std::max(abs(rawSidechainSignal[i]), 1e-6);
                        //rawSidechainSignal[i] = abs(rawSidechainSignal[i]);
                        //rawSidechainSignal[i] = smooth::ABS(rawSidechainSignal[i], 0.0001);

                        //outputs[ch][i] = (SampleType)rawSidechainSignal[i];
                        //outputs[ch][i] = std::abs(inputs[ch][i]);
                    }
                }
                if (false) {
                    for (int i = 0; i < sampleFrames; i++) {
                        // two paths for 0 and +90
                        for (int p = 0; p < 2; p++) {
                            double ret1 = c[p][0] * (rawSidechainSignal[i] + p1_sy_1[p]) - p1_sx_1[p];
                            p1_sx_1[p] = p1_sx_2[p];
                            p1_sx_2[p] = rawSidechainSignal[i];
                            p1_sy_1[p] = p1_sy_2[p];
                            p1_sy_2[p] = ret1;
                            double ret2 = c[p][1] * (ret1 + p2_sy_1[p]) - p2_sx_1[p];
                            p2_sx_1[p] = p2_sx_2[p];
                            p2_sx_2[p] = ret1; // inputs[ch][i];
                            p2_sy_1[p] = p2_sy_2[p];
                            p2_sy_2[p] = ret2;
                            double ret3 = c[p][2] * (ret2 + p3_sy_1[p]) - p3_sx_1[p];
                            p3_sx_1[p] = p3_sx_2[p];
                            p3_sx_2[p] = ret2; // inputs[ch][i];
                            p3_sy_1[p] = p3_sy_2[p];
                            p3_sy_2[p] = ret3;
                            double ret4 = c[p][3] * (ret3 + p4_sy_1[p]) - p4_sx_1[p];
                            p4_sx_1[p] = p4_sx_2[p];
                            p4_sx_2[p] = ret3; // inputs[ch][i];
                            p4_sy_1[p] = p4_sy_2[p];
                            p4_sy_2[p] = ret4;
                        }
                        double env = sqrt(p4_sy_1[1] * p4_sy_1[1] + p4_sy_2[0] * p4_sy_2[0]); // *0.637?
                        rawSidechainSignal[i] = env;
                    }
                }
                //return;

                // upsampling
                for (int ch = 0; ch < 1; ch++) {
                    for (int i = 0; i < sampleFrames; i++) {
                        // OS - upsample
                        up[ch].acc();
                        up[ch].buff[up[ch].now] = (double)rawSidechainSignal[i];
                        //up[ch].buff[up[ch].now] = (double)inputs[ch][i];
                        double in0 = 0.0;
                        for (int i = 0; i < tap_hm; i++) {
                            double a = up[ch].buff[up[ch].get_nth(i)];
                            double b = up[ch].buff[up[ch].get_nth(fir_size - 1 - i)];
                            in0 += up[ch].coef[i] * (a + b);
                        }   in0 += up[ch].coef[tap_hm] * up[ch].buff[up[ch].get_nth(tap_hm)];
                        up[ch].acc();
                        up[ch].buff[up[ch].now] = 0.0;
                        double in1 = 0.0;
                        for (int i = 0; i < tap_hm; i++) {
                            double a = up[ch].buff[up[ch].get_nth(i)];
                            double b = up[ch].buff[up[ch].get_nth(fir_size - 1 - i)];
                            in1 += up[ch].coef[i] * (a + b);
                        }   in1 += up[ch].coef[tap_hm] * up[ch].buff[up[ch].get_nth(tap_hm)];

                        upsample_buff[ch][(2 * i) + 0] = in0;
                        upsample_buff[ch][(2 * i) + 1] = in1;
                    }
                }

                //int ch = 0;
                /*
                // Get max l/r amplitude values and fill sidechain signal
                for (int ch = 0; ch < numChannels; ch++) {
                    for (int i = 0; i < 2*sampleFrames; i++) {
                        upsample_buff[0][i] = (abs(upsample_buff[0][i]) < abs(upsample_buff[ch][i]))
                            ? upsample_buff[0][i] : upsample_buff[ch][i];
                    }
                }
                */
                for (int i = 0; i < 2 * sampleFrames; i++)
                    upsample_buff[0][i] = smooth::ABS(upsample_buff[0][i], 0.0001);
                ballistics.processCrestFactor(upsample_buff[0].data(), 2 * sampleFrames);
                gainComputer.applyCompressionToBuffer(upsample_buff[0].data(), 2 * sampleFrames);
                ballistics.applyBallistics(upsample_buff[0].data(), 2 * sampleFrames);

                // Do lookahead if enabled
                if (true)//LookAheadEnable
                {
                    // Delay input buffer
                    //delay.process(buffer);
                    int ch = 0;
                    int latency = 128; // (int)lookaheadDelay;
                    if (latency != delayline[ch].size()) {
                        int32 diff = latency - (int32)delayline[ch].size();
                        if (diff > 0)
                            for (int i = 0; i < diff; i++) delayline[ch].push_back(0.0);
                        else
                            for (int i = 0; i < -diff; i++) delayline[ch].pop_front();
                    }

                    for (int i = 0; i < 2 * sampleFrames; i++) {
                        double ahead = upsample_buff[0][i];
                        double delayed = delayline[ch].front();
                        delayline[ch].pop_front();
                        delayline[ch].push_back(upsample_buff[0][i]);

                        double acc = std::accumulate(delayline[ch]->begin(), delayline[ch]->end(), 0.0);
                        acc /= delayline[ch].size();
                        // upsample_buff[0][i] = acc;

                        if (ahead > acc) upsample_buff[0][i] = acc;
                        else upsample_buff[0][i] = acc;

                        //if (acc < delayed) upsample_buff[0][i] = acc;
                        //else upsample_buff[0][i] = delayed;
                    }


                    // Gain reduction must be smoothed with side chain window width
                    // Maybe, crossfading with LAH / ORIG only when LAH > ORIG

                    // Process side-chain (delay + gain reduction fade in)

                }

                for (int ch = 0; ch < 1; ch++) {
                    for (int i = 0; i < sampleFrames; i++) {
                        // OS - downsample
                        //memmove(dn[ch].buff + 1, dn[ch].buff, sizeof(double) * (maxTap)-1);
                        dn[ch].acc();
                        dn[ch].buff[dn[ch].now] = upsample_buff[ch][(2 * i) + 0];

                        double out = 0.0;

                        for (int i = 0; i < tap_hm; i++) {
                            double a = dn[ch].buff[dn[ch].get_nth(i)];
                            double b = dn[ch].buff[dn[ch].get_nth(fir_size - 1 - i)];
                            out += dn[ch].coef[i] * (a + b);
                        }   out += dn[ch].coef[tap_hm] * dn[ch].buff[dn[ch].get_nth(tap_hm)];

                        //memmove(dn[ch].buff + 1, dn[ch].buff, sizeof(double) * (maxTap)-1);
                        dn[ch].acc();
                        dn[ch].buff[dn[ch].now] = upsample_buff[ch][(2 * i) + 1];

                        rawSidechainSignal[i] = (SampleType)out;
                    }
                }

                /* - original

                // Calculate crest factor on max. amplitude values of input buffer
                ballistics.processCrestFactor(rawSidechainSignal, sampleFrames);

                // Compute attenuation - converts side-chain signal from linear to logarithmic domain
                gainComputer.applyCompressionToBuffer(rawSidechainSignal, sampleFrames);
                // Now rawSidechainSignal has ideal GR dB.
                // Smooth attenuation - still logarithmic
                ballistics.applyBallistics(rawSidechainSignal, sampleFrames);

                */

                // Get minimum = max. gain reduction from side chain buffer
                double min = 0;
                for (int i = 0; i < sampleFrames; i++)
                    if (rawSidechainSignal[i] < min)
                        min = rawSidechainSignal[i];
                maxGainReduction = min;



                // Calculate auto makeup
                autoMakeup = calculateAutoMakeup(rawSidechainSignal, sampleFrames);

                // Add makeup gain and convert side-chain to linear domain
                for (int i = 0; i < sampleFrames; i++)
                    sidechainSignal[i] = Decibels::decibelsToGain(sidechainSignal[i] + makeup + autoMakeup);



                for (int ch = 0; ch < numChannels; ch++)
                {
                    int latency = 63 + 64; // (int)lookaheadDelay;
                    if (latency != latency_q[ch].size()) {
                        int32 diff = latency - (int32)latency_q[ch].size();
                        if (diff > 0) for (int i = 0; i < diff; i++) latency_q[ch].push(0.0);
                        else for (int i = 0; i < -diff; i++) latency_q[ch].pop();
                    }

                    for (int j = 0; j < sampleFrames; j++) {
                        // Copy buffer to original signal
                        // originalSignal[ch][j] = (double)inputs[ch][j];

                        originalSignal[ch][j] = latency_q[ch].front();
                        latency_q[ch].pop();
                        latency_q[ch].push((double)inputs[ch][j]);


                        // Multiply attenuation with buffer - apply compression
                        double t = originalSignal[ch][j] * rawSidechainSignal[j];

                        // Mix dry & wet signal
                        outputs[ch][j] = (SampleType)(t * mix + originalSignal[ch][j] * (1.0 - mix));
                    }
                }
            }
        };

    private:
        template <typename SampleType>
        inline void applyInputGain(SampleType* buffer, int numSamples)
        {
            for (int i = 0; i < numSamples; i++) buffer[i] *= Decibels::decibelsToGain(input);
            /*
            if (prevInput == input)
                buffer.applyGain(0, numSamples, Decibels::decibelsToGain(prevInput));
            else
            {
                buffer.applyGainRamp(0, numSamples, Decibels::decibelsToGain(prevInput), Decibels::decibelsToGain(input));
                prevInput = input;
            }
            */
        };
        inline double calculateAutoMakeup(const double* src, int numSamples)
        {
            double sum = 0.0;
            for (int i = 0; i < numSamples; i++) {
                sum += src[i];
            }

            smoothedAutoMakeup.process(-sum / static_cast<double>(numSamples));
            return autoMakeupEnabled ? static_cast<double>(smoothedAutoMakeup.getState()) : 0.0f;
        };

        //Directly initialize process spec to avoid debugging problems
        double _sampleRate;
        double _maximumBlockSize;

        std::vector<double> originalSignal[2];
        std::vector<double> sidechainSignal;
        double* rawSidechainSignal{ nullptr };

        LevelDetector ballistics;
        GainComputer  gainComputer;
        SmoothingFilter smoothedAutoMakeup;

        double lookaheadDelay{ 0.005 };
        double input{ 0.0 };
        double prevInput{ 0.0 };
        double makeup{ 0.0 };
        double autoMakeup{ 0.0 };
        bool   bypassed{ false };
        bool   autoMakeupEnabled{ false };
        bool   LookAheadEnable{ true };
        double mix{ 1.0 };
        double maxGainReduction{ 0.0 };

        Flt up[2];
        Flt dn[2];
        std::vector<double> upsample_buff[2];
        std::vector<double> dnsample_buff[2];
        std::deque<double> delayline[2];
        const int fir_size = 125;
        const int tap_hm = (fir_size - 1) / 2;
        std::queue<double> latency_q[2];

        double c[2][4] = {
        {0.16514909355907719801,0.73982901254452670958,0.94794090632917971107,0.99120971270525837227},
        {0.48660436861367767358,0.88077943527246449484,0.97793125561632343601,0.99767386185073303473}
        };
        double p1_sx_1[2] = { 0.0, };
        double p1_sx_2[2] = { 0.0, };
        double p1_sy_1[2] = { 0.0, };
        double p1_sy_2[2] = { 0.0, };
        double p2_sx_1[2] = { 0.0, };
        double p2_sx_2[2] = { 0.0, };
        double p2_sy_1[2] = { 0.0, };
        double p2_sy_2[2] = { 0.0, };
        double p3_sx_1[2] = { 0.0, };
        double p3_sx_2[2] = { 0.0, };
        double p3_sy_1[2] = { 0.0, };
        double p3_sy_2[2] = { 0.0, };
        double p4_sx_1[2] = { 0.0, };
        double p4_sx_2[2] = { 0.0, };
        double p4_sy_1[2] = { 0.0, };
        double p4_sy_2[2] = { 0.0, };
    };


    /* Basic envelope-follwer, to track peak & rms signal level with configurable decay time*/
    class LevelEnvelopeFollower
    {
    public:
        LevelEnvelopeFollower() = default;

        // Prepares envelope follower with given sample rate and recalculates decayInSamples
        // as well as the peak/rms coefficient
        void prepare(const double& fs)
        {
            sampleRate = fs;

            peakDecayInSamples = static_cast<int>(peakDecayInSeconds * sampleRate);
            peakDecay = 1.0 - 1.0 / static_cast<double>(peakDecayInSamples);

            rmsDecayInSamples = static_cast<int>(rmsDecayInSeconds * sampleRate);
            rmsDecay = 1.0 - 1.0 / static_cast<double>(rmsDecayInSamples);

            double attackTimeInSeconds = 0.0; //Time it takes to reach 1-1/e = 0.63
            alphaAttack = 0.0;// exp(-1.0 / (sampleRate * attackTimeInSeconds)); //aA = e^(-1/TA*fs)

            double releaseTimeInSeconds = peakDecayInSeconds; //Time it takes to reach 1 - (1-1/e) = 0.37
            alphaRelease = exp(-1.0 / (sampleRate * releaseTimeInSeconds)); //aR = e^(-1/TR*fs)

            currMaxPeak[0] = Decibels::gainToDecibels(0.0);
            currMaxPeak[1] = Decibels::gainToDecibels(0.0);
            currMaxRMS[0] = Decibels::gainToDecibels(0.0);
            currMaxRMS[1] = Decibels::gainToDecibels(0.0);
            state[0] = Decibels::gainToDecibels(0.0);
            state[1] = Decibels::gainToDecibels(0.0);
        };

        // Set peak decay
        void setPeakDecay(double dc)
        {
            peakDecayInSeconds = dc;
            prepare(sampleRate);
        };

        // Set rms decay
        void setRmsDecay(double dc) {
            rmsDecayInSeconds = dc;
            prepare(sampleRate);
        };

        // Updates peak envelope follower from given audio buffer
        template <typename SampleType>
        void updatePeak(const SampleType* const* channelData, int numChannels, int numSamples)
        {
            if (numChannels > 0 && numSamples > 0) {
                for (int ch = 0; ch < numChannels; ch++) {
                    for (int i = 0; i < numSamples; i++) {
                        double in = Decibels::gainToDecibels(std::abs(channelData[ch][i]));

                        if (in > state[ch])
                            state[ch] = alphaAttack * state[ch] + (1 - alphaAttack) * in;
                        else
                            state[ch] = alphaRelease * state[ch] + (1 - alphaRelease) * in;

                        currMaxPeak[ch] = (state[ch]); //y_L
                    }
                    currMaxPeak[ch] = std::min(6.0, (currMaxPeak[ch]));
                }
            }
        };

        // Updates rms envelope follower from given audio buffer
        void updateRMS(const double* const* channelData, int numChannels, int numSamples)
        {
            if (numChannels > 0 && numSamples > 0) {
                for (int ch = 0; ch < numChannels; ch++) {
                    for (int i = 0; i < numSamples; i++) {
                        double sum = std::abs(channelData[ch][i]);
                        sum *= sum;

                        if (sum > currMaxRMS[ch])
                            currMaxRMS[ch] = sum * rmsDecay;
                        else if (currMaxRMS[ch] > 0.001f)
                            currMaxRMS[ch] *= peakDecay;
                        else currMaxRMS[ch] = 0.0f;
                    }
                }
            }
        };

        // Gets current peak, call after updatePeak
        double getPeak(int channel) {
            return currMaxPeak[channel];
        };

        // Gets current rms, vall after updateRMS
        double getRMS(int channel) {
            return sqrt(currMaxRMS[channel]);
        };

    private:
        double currMaxPeak[2] = { 0.0, 0.0 };
        double currMaxRMS[2] = { 0.0, 0.0 };
        double peakDecay{ 0.99992 };
        double rmsDecay{ 0.95 };
        double peakDecayInSeconds{ 0.5 };
        double rmsDecayInSeconds{ 0.0 };

        int peakDecayInSamples{ 0 };
        int rmsDecayInSamples{ 0 };

        double state[2] = { 0.0, 0.0 };
        double alphaAttack = 0.0;
        double alphaRelease = 0.0;
        double sampleRate{ 0.0 };
    };


    using namespace Steinberg;

//------------------------------------------------------------------------
static constexpr Vst::DataExchangeBlock 
    InvalidDataExchangeBlock = { nullptr, 0, Vst::InvalidDataExchangeBlockID };

//------------------------------------------------------------------------
//  comp_Processor
//------------------------------------------------------------------------
class comp_Processor : public Vst::AudioEffect
{
public:
	comp_Processor ();
	~comp_Processor () SMTG_OVERRIDE;

    // Create function
	static FUnknown* createInstance (void* /*context*/) 
	{ 
		return (Steinberg::Vst::IAudioProcessor*)new comp_Processor; 
	}

	//--- ---------------------------------------------------------------------
	// AudioEffect overrides:
	//--- ---------------------------------------------------------------------
	/** Called at first after constructor */
	Steinberg::tresult PLUGIN_API initialize (Steinberg::FUnknown* context) SMTG_OVERRIDE;
	
	/** Called at the end before destructor */
	Steinberg::tresult PLUGIN_API terminate () SMTG_OVERRIDE;

    /** Try to set (host => plug-in) a wanted arrangement for inputs and outputs. */
	Steinberg::tresult PLUGIN_API setBusArrangements(
		Steinberg::Vst::SpeakerArrangement* inputs, Steinberg::int32 numIns,
		Steinberg::Vst::SpeakerArrangement* outputs, Steinberg::int32 numOuts
	) SMTG_OVERRIDE;

    /** Connects this instance with another connection point. */
    Steinberg::tresult PLUGIN_API connect(Steinberg::Vst::IConnectionPoint* other) override;

    /** Disconnects a given connection point from this. */
    Steinberg::tresult PLUGIN_API disconnect(Steinberg::Vst::IConnectionPoint* other) override;
	
	/** Switch the Plug-in on/off */
	Steinberg::tresult PLUGIN_API setActive (Steinberg::TBool state) SMTG_OVERRIDE;

	/** Will be called before any process call */
	Steinberg::tresult PLUGIN_API setupProcessing (Steinberg::Vst::ProcessSetup& newSetup) SMTG_OVERRIDE;
	
	/** Asks if a given sample size is supported see SymbolicSampleSizes. */
	Steinberg::tresult PLUGIN_API canProcessSampleSize (Steinberg::int32 symbolicSampleSize) SMTG_OVERRIDE;

    /** Gets the current Latency in samples. */
    Steinberg::uint32  PLUGIN_API getLatencySamples() SMTG_OVERRIDE;

	/** Here we go...the process call */
	Steinberg::tresult PLUGIN_API process (Steinberg::Vst::ProcessData& data) SMTG_OVERRIDE;
		
	/** For persistence */
	Steinberg::tresult PLUGIN_API setState (Steinberg::IBStream* state) SMTG_OVERRIDE;
	Steinberg::tresult PLUGIN_API getState (Steinberg::IBStream* state) SMTG_OVERRIDE;

    double gainReduction = 0.0;
    double currentInput[2] = { -std::numeric_limits<double>::infinity(), };
    double currentOutput[2] = { -std::numeric_limits<double>::infinity(), };

    template <typename SampleType>
    void processComp(
        SampleType** inputs,
        SampleType** outputs,
        Steinberg::int32 numChannels,
        Steinberg::Vst::SampleRate getSampleRate,
        Steinberg::int32 sampleFrames
    );

//------------------------------------------------------------------------
protected:

    double Attack_LOG_MAX = log(50.0 / 0.1);
    double Release_LOG_MAX = log(5000.0 / 5.0);

	bool            bBypass          = false;
	Vst::ParamValue fZoom            = 2.0 / 6.0;
    Vst::ParamValue fThreshold       = (-12.0 - (-60.0)) / (0.0 - (-60.0));
    Vst::ParamValue fRatio           = (4.0 - 1.0) / (20.0 - 1.0);
    Vst::ParamValue fKnee            = (5.0 - 0.0) / (20.0 - 0.0);
    Vst::ParamValue fAttack          = log(1.0 / 0.1) / Attack_LOG_MAX;
	Vst::ParamValue fRelease         = log(150.0 / 5.0) / Release_LOG_MAX;
    Vst::ParamValue fLookAhead       = (0.0 - 0.0) / (512.0 - 0.0);
    bool            bLookAheadEnable = false;
	Vst::ParamValue fMakeup          = (0.0 - (-12.0)) / (12.0 - (-12.0));
	Vst::ParamValue fMix             = 1.0;

    Compressor compressor;
    LevelEnvelopeFollower inLevelFollower;
    LevelEnvelopeFollower outLevelFollower;


    Flt up[2];
    Flt dn[2];
    const int fir_size = 125;
    const int tap_hm = (fir_size - 1) / 2;
    std::vector<double> tt[2];
    std::vector<double> dd[2];
    double* ttt[2]{ nullptr, };
    double* ddd[2]{ nullptr, };

    void acquireNewExchangeBlock();
    std::unique_ptr<Steinberg::Vst::DataExchangeHandler> dataExchange;
    Steinberg::Vst::DataExchangeBlock currentExchangeBlock{ InvalidDataExchangeBlock };
    uint16_t numChannels{ 0 };

    int32 latency = 63;

    double coeff = exp(-1.0 / (0.1 * 0.001 * 48000.0));
    double icoef = 1.0 - coeff;
    double rms_s = 0.0, rms_c = 0.0, rms = 0.0;
    double peak = 0.0, pp = 0.0;

    // https://yehar.com/blog/?p=368
    // https://dsp.stackexchange.com/questions/37411/iir-hilbert-transformer
    // https://forum.juce.com/t/90-phase-shift/6822/5
    double c[2][4] = {
        {0.16514909355907719801,0.73982901254452670958,0.94794090632917971107,0.99120971270525837227},
        {0.48660436861367767358,0.88077943527246449484,0.97793125561632343601,0.99767386185073303473}
    };
    double p1_sx_1[2] = { 0.0, };
    double p1_sx_2[2] = { 0.0, };
    double p1_sy_1[2] = { 0.0, };
    double p1_sy_2[2] = { 0.0, };
    double p2_sx_1[2] = { 0.0, };
    double p2_sx_2[2] = { 0.0, };
    double p2_sy_1[2] = { 0.0, };
    double p2_sy_2[2] = { 0.0, };
    double p3_sx_1[2] = { 0.0, };
    double p3_sx_2[2] = { 0.0, };
    double p3_sy_1[2] = { 0.0, };
    double p3_sy_2[2] = { 0.0, };
    double p4_sx_1[2] = { 0.0, };
    double p4_sx_2[2] = { 0.0, };
    double p4_sy_1[2] = { 0.0, };
    double p4_sy_2[2] = { 0.0, };

    double sx_1[2] = { 0.0, };
    double sx_2[2] = { 0.0, };
    double sy_1[2] = { 0.0, };
    double sy_2[2] = { 0.0, };

    std::queue<double> qq[2];
};

//------------------------------------------------------------------------
} // namespace yg331
