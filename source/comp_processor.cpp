//------------------------------------------------------------------------
// Copyright(c) 2024 yg331.
//------------------------------------------------------------------------

#include "comp_processor.h"
//#include "comp_cids.h"

#include "base/source/fstreamer.h"
#include "pluginterfaces/vst/ivstparameterchanges.h"
#include "public.sdk/source/vst/vstaudioprocessoralgo.h"
#include "public.sdk/source/vst/vsthelpers.h"

using namespace Steinberg;

namespace yg331 {
//------------------------------------------------------------------------
// comp_Processor
//------------------------------------------------------------------------
comp_Processor::comp_Processor ()
{
	//--- set the wanted controller for our processor
	setControllerClass (kcomp_ControllerUID);
}

//------------------------------------------------------------------------
comp_Processor::~comp_Processor ()
{}

//------------------------------------------------------------------------
tresult PLUGIN_API comp_Processor::initialize (FUnknown* context)
{
	// Here the Plug-in will be instantiated
	
	//---always initialize the parent-------
	tresult result = AudioEffect::initialize (context);
	// if everything Ok, continue
	if (result != kResultOk)
	{
		return result;
	}

	//--- create Audio IO ------
	addAudioInput (STR16 ("Stereo In"), Steinberg::Vst::SpeakerArr::kStereo);
	addAudioOutput (STR16 ("Stereo Out"), Steinberg::Vst::SpeakerArr::kStereo);

	/* If you don't need an event bus, you can remove the next line */
	//saddEventInput (STR16 ("Event In"), 1);
	
	for (int ch = 0; ch < 2; ch++) {
		up[ch].size = fir_size;
		dn[ch].size = fir_size;
		Kaiser::calcFilter(96000.0, 0.0, 24000.0, fir_size, 200.0, up[ch].coef);
		Kaiser::calcFilter(96000.0, 0.0, 24000.0, fir_size, 200.0, dn[ch].coef);

		for (int i = 0; i < fir_size; i++) dn[ch].coef[i] *= 2.0;
	}

	return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API comp_Processor::terminate ()
{
	// Here the Plug-in will be de-instantiated, last possibility to remove some memory!
	
	//---do not forget to call parent ------
	return AudioEffect::terminate ();
}

//------------------------------------------------------------------------
tresult PLUGIN_API comp_Processor::setBusArrangements(
	Steinberg::Vst::SpeakerArrangement* inputs, Steinberg::int32 numIns, 
	Steinberg::Vst::SpeakerArrangement* outputs, Steinberg::int32 numOuts)
{
	if (numIns == 1 && numOuts == 1)
	{
		// the host wants Mono => Mono (or 1 channel -> 1 channel)
		if (Vst::SpeakerArr::getChannelCount(inputs[0]) == 1 &&
			Vst::SpeakerArr::getChannelCount(outputs[0]) == 1)
		{
			auto* bus = FCast<Vst::AudioBus>(audioInputs.at(0));
			if (bus)
			{
				// check if we are Mono => Mono, if not we need to recreate the busses
				if (bus->getArrangement() != inputs[0])
				{
					getAudioInput(0)->setArrangement(inputs[0]);
					getAudioInput(0)->setName(STR16("Mono In"));
					getAudioOutput(0)->setArrangement(outputs[0]);
					getAudioOutput(0)->setName(STR16("Mono Out"));
				}
				return kResultOk;
			}
		}
		// the host wants something else than Mono => Mono,
		// in this case we are always Stereo => Stereo
		else
		{
			auto* bus = FCast<Vst::AudioBus>(audioInputs.at(0));
			if (bus)
			{
				tresult result = kResultFalse;

				// the host wants 2->2 (could be LsRs -> LsRs)
				if (Vst::SpeakerArr::getChannelCount(inputs[0]) == 2 &&
					Vst::SpeakerArr::getChannelCount(outputs[0]) == 2)
				{
					getAudioInput(0)->setArrangement(inputs[0]);
					getAudioInput(0)->setName(STR16("Stereo In"));
					getAudioOutput(0)->setArrangement(outputs[0]);
					getAudioOutput(0)->setName(STR16("Stereo Out"));
					result = kResultTrue;
				}
				// the host want something different than 1->1 or 2->2 : in this case we want stereo
				else if (bus->getArrangement() != Vst::SpeakerArr::kStereo)
				{
					getAudioInput(0)->setArrangement(Vst::SpeakerArr::kStereo);
					getAudioInput(0)->setName(STR16("Stereo In"));
					getAudioOutput(0)->setArrangement(Vst::SpeakerArr::kStereo);
					getAudioOutput(0)->setName(STR16("Stereo Out"));
					result = kResultFalse;
				}

				return result;
			}
		}
	}
	return kResultFalse;
}


//------------------------------------------------------------------------
tresult PLUGIN_API comp_Processor::connect(Vst::IConnectionPoint* other)
{
	auto result = Vst::AudioEffect::connect (other);
	if (result == kResultTrue)
	{
		auto configCallback = [this] (
			Vst::DataExchangeHandler::Config& config,
		    const Vst::ProcessSetup& setup
		) 
		{
			Vst::SpeakerArrangement arr;
			getBusArrangement (Vst::BusDirections::kInput, 0, arr);
			numChannels = static_cast<uint16_t> (Vst::SpeakerArr::getChannelCount (arr));
			auto sampleSize = sizeof (float);

			config.blockSize = numChannels * sampleSize + sizeof (DataBlock);
			config.numBlocks = 2;
			config.alignment = 32;
			config.userContextID = 0;
			return true;
		};

		dataExchange = std::make_unique<Vst::DataExchangeHandler> (this, configCallback);
		dataExchange->onConnect (other, getHostContext ());
	}
	return result;
}

//------------------------------------------------------------------------
tresult PLUGIN_API comp_Processor::disconnect(Vst::IConnectionPoint* other)
{
	if (dataExchange)
	{
		dataExchange->onDisconnect(other);
		dataExchange.reset();
	}
	return AudioEffect::disconnect(other);
}


//------------------------------------------------------------------------
tresult PLUGIN_API comp_Processor::setActive (TBool state)
{
	//--- called when the Plug-in is enable/disable (On/Off) -----
	
	if (state)
		dataExchange->onActivate(processSetup);
	else
		dataExchange->onDeactivate();
	
	return AudioEffect::setActive (state);
}

//------------------------------------------------------------------------
void comp_Processor::acquireNewExchangeBlock()
{
	currentExchangeBlock = dataExchange->getCurrentOrNewBlock();
	if (auto block = toDataBlock(currentExchangeBlock))
	{
		block->sampleRate = static_cast<uint32_t> (processSetup.sampleRate);
		//block->numChannels = numChannels;
		//block->sampleSize = sizeof(float);
		block->numSamples = 0;
		block->inL = currentInput[0];
		block->inR = currentInput[1];
		block->outL = currentOutput[0];
		block->outR = currentOutput[1];
		block->gR = gainReduction;
	}
}

uint32 PLUGIN_API comp_Processor::getLatencySamples()
{
	return latency; // 63 + 128;
}

//------------------------------------------------------------------------
tresult PLUGIN_API comp_Processor::process (Vst::ProcessData& data)
{
	Vst::IParameterChanges* paramChanges = data.inputParameterChanges;

	if (paramChanges)
	{
		int32 numParamsChanged = paramChanges->getParameterCount();

		for (int32 index = 0; index < numParamsChanged; index++)
		{
			Vst::IParamValueQueue* paramQueue = paramChanges->getParameterData(index);

			if (paramQueue)
			{
				Vst::ParamValue value;
				int32 sampleOffset;
				int32 numPoints = paramQueue->getPointCount();

				/*/*/
				if (paramQueue->getPoint(numPoints - 1, sampleOffset, value) == kResultTrue) {
					switch (paramQueue->getParameterId()) {
					case kParamBypass:          bBypass    = (value > 0.5f); break;
					case kParamZoom:            fZoom      = value; break;
					case kParamThreshold:       fThreshold = value; break;
					case kParamRatio:           fRatio     = value; break;
					case kParamKnee:            fKnee      = value; break;
					case kParamAttack:          fAttack    = value; break;
					case kParamRelease:         fRelease   = value; break;
					case kParamLookAhead:       fLookAhead = value; break;
					case kParamLookAheadEnable: bLookAheadEnable = (value > 0.5f); break;
					case kParamMakeup:          fMakeup    = value; break;
					case kParamMix:             fMix       = value; break;
					}
				}
			}
		}
	}

	if (data.numInputs == 0 || data.numOutputs == 0)
	{
		// nothing to do
		return kResultOk;
	}

	// (simplification) we suppose in this example that we have the same input channel count than
	// the output
	int32 numChannels = data.inputs[0].numChannels;

	//---get audio buffers----------------
	uint32 sampleFramesSize = getSampleFramesSizeInBytes(processSetup, data.numSamples);
	void** in = getChannelBuffersPointer(processSetup, data.inputs[0]);
	void** out = getChannelBuffersPointer(processSetup, data.outputs[0]);
	Vst::SampleRate getSampleRate = processSetup.sampleRate;

	//---check if silence---------------
	// check if all channel are silent then process silent
	if (data.inputs[0].silenceFlags == Vst::getChannelMask(data.inputs[0].numChannels))
	{
		// mark output silence too (it will help the host to propagate the silence)
		data.outputs[0].silenceFlags = data.inputs[0].silenceFlags;

		// the plug-in has to be sure that if it sets the flags silence that the output buffer are
		// clear
		for (int32 i = 0; i < numChannels; i++)
		{
			// do not need to be cleared if the buffers are the same (in this case input buffer are
			// already cleared by the host)
			if (in[i] != out[i])
			{
				memset(out[i], 0, sampleFramesSize);
			}
		}
	}
	else {

		data.outputs[0].silenceFlags = data.inputs[0].silenceFlags;

		//---in bypass mode outputs should be like inputs-----
		if (bBypass)
		{
			if (data.symbolicSampleSize == Vst::kSample32) {
				//latencyBypass<Vst::Sample32>((Vst::Sample32**)in, (Vst::Sample32**)out, numChannels, getSampleRate, data.numSamples);
			}
			else if (data.symbolicSampleSize == Vst::kSample64) {
				//latencyBypass<Vst::Sample64>((Vst::Sample64**)in, (Vst::Sample64**)out, numChannels, getSampleRate, data.numSamples);
			}
		}
		else {
			if (data.symbolicSampleSize == Vst::kSample32) {
				processComp<Vst::Sample32>((Vst::Sample32**)in, (Vst::Sample32**)out, numChannels, getSampleRate, data.numSamples);
			}
			else {
				processComp<Vst::Sample64>((Vst::Sample64**)in, (Vst::Sample64**)out, numChannels, getSampleRate, data.numSamples);
			}
		}
	}

	return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API comp_Processor::setupProcessing (Vst::ProcessSetup& newSetup)
{
	//--- called before any processing ----
	compressor.prepare(newSetup.sampleRate, newSetup.maxSamplesPerBlock);
	inLevelFollower.prepare(newSetup.sampleRate);
	outLevelFollower.prepare(newSetup.sampleRate);
	inLevelFollower.setPeakDecay(3.0);
	outLevelFollower.setPeakDecay(3.0);
	tt[0].resize(2 * newSetup.maxSamplesPerBlock, 0.0);
	tt[1].resize(2 * newSetup.maxSamplesPerBlock, 0.0);
	dd[0].resize(2 * newSetup.maxSamplesPerBlock, 0.0);
	dd[1].resize(2 * newSetup.maxSamplesPerBlock, 0.0);
	ttt[0] = tt[0].data();
	ttt[1] = tt[1].data();
	ddd[0] = dd[0].data();
	ddd[1] = dd[1].data();
	return AudioEffect::setupProcessing (newSetup);
}

//------------------------------------------------------------------------
tresult PLUGIN_API comp_Processor::canProcessSampleSize (int32 symbolicSampleSize)
{
	// by default kSample32 is supported
	if (symbolicSampleSize == Vst::kSample32)
		return kResultTrue;

	// disable the following comment if your processing support kSample64
	if (symbolicSampleSize == Vst::kSample64)
		return kResultTrue;

	return kResultFalse;
}

//------------------------------------------------------------------------
tresult PLUGIN_API comp_Processor::setState (IBStream* state)
{
	// called when we load a preset, the model has to be reloaded
	IBStreamer streamer (state, kLittleEndian);

	int32           savedBypass          = 0;
	Vst::ParamValue savedZoom            = 0.0;

	Vst::ParamValue savedThreshold       = 0.0;
	Vst::ParamValue savedRatio           = 0.0;
	Vst::ParamValue savedKnee            = 0.0;
	Vst::ParamValue savedAttack          = 0.0;
	Vst::ParamValue savedRelease         = 0.0;
	Vst::ParamValue savedLookAhead       = 0.0;
	int32           savedLookAheadEnable = 0.0;

	Vst::ParamValue savedMakeup          = 0.0;
	Vst::ParamValue savedMix             = 0.0;

	if (streamer.readInt32 (savedBypass         ) == false) return kResultFalse;
	if (streamer.readDouble(savedZoom           ) == false) return kResultFalse;
	if (streamer.readDouble(savedThreshold      ) == false) return kResultFalse;
	if (streamer.readDouble(savedRatio          ) == false) return kResultFalse;
	if (streamer.readDouble(savedKnee           ) == false) return kResultFalse;
	if (streamer.readDouble(savedAttack         ) == false) return kResultFalse;
	if (streamer.readDouble(savedRelease        ) == false) return kResultFalse;
	if (streamer.readDouble(savedLookAhead      ) == false) return kResultFalse;
	if (streamer.readInt32 (savedLookAheadEnable) == false) return kResultFalse;
	if (streamer.readDouble(savedMakeup         ) == false) return kResultFalse;
	if (streamer.readDouble(savedMix            ) == false) return kResultFalse;

	bBypass          = savedBypass > 0;
	fZoom            = savedZoom;
	fThreshold       = savedThreshold;
	fRatio           = savedRatio;
	fKnee            = savedKnee;
	fAttack          = savedAttack;
	fRelease         = savedRelease;
	fLookAhead       = savedLookAhead;
	bLookAheadEnable = savedLookAheadEnable;
	fMakeup          = savedMakeup;
	fMix             = savedMix;

	if (Vst::Helpers::isProjectState(state) == kResultTrue)
	{
		// we are in project loading context...

		// Example of using the IStreamAttributes interface
		FUnknownPtr<Vst::IStreamAttributes> stream(state);
		if (stream)
		{
			if (Vst::IAttributeList* list = stream->getAttributes())
			{
				// get the full file path of this state
				Vst::TChar fullPath[1024];
				memset(fullPath, 0, 1024 * sizeof(Vst::TChar));
				if (list->getString(Vst::PresetAttributes::kFilePathStringType, fullPath,
					1024 * sizeof(Vst::TChar)) == kResultTrue)
				{
					// here we have the full path ...
				}
			}
		}
	}
	return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API comp_Processor::getState (IBStream* state)
{
	// here we need to save the model
	IBStreamer streamer (state, kLittleEndian);

	streamer.writeInt32 (bBypass ? 1 : 0);
	streamer.writeDouble(fZoom);
	streamer.writeDouble(fThreshold);
	streamer.writeDouble(fRatio);
	streamer.writeDouble(fKnee);
	streamer.writeDouble(fAttack);
	streamer.writeDouble(fRelease);
	streamer.writeDouble(fLookAhead);
	streamer.writeInt32 (bLookAheadEnable ? 1 : 0);
	streamer.writeDouble(fMakeup);
	streamer.writeDouble(fMix);

	return kResultOk;
}


template <typename SampleType>
void comp_Processor::processComp(
	SampleType** inputs,
	SampleType** outputs,
	Steinberg::int32 numChannels,
	Steinberg::Vst::SampleRate getSampleRate,
	Steinberg::int32 sampleFrames)
{
	if (false) {
		for (int ch = 0; ch < numChannels; ch++) {
			int latency = 1;
			if (latency != qq[ch].size()) {
				int32 diff = latency - (int32)qq[ch].size();
				if (diff > 0) for (int i = 0; i < diff; i++) qq[ch].push(0.0);
				else for (int i = 0; i < -diff; i++) qq[ch].pop();
			}
			for (int i = 0; i < sampleFrames; i++) {
				qq[ch].push(inputs[ch][i]);
				double t = inputs[ch][i];
				double ret1 = c[ch][0] * (t + p1_sy_1[ch]) - p1_sx_1[ch];
				p1_sx_1[ch] = p1_sx_2[ch];
				p1_sx_2[ch] = t;
				p1_sy_1[ch] = p1_sy_2[ch];
				p1_sy_2[ch] = ret1;
				double ret2 = c[ch][1] * (ret1 + p2_sy_1[ch]) - p2_sx_1[ch];
				p2_sx_1[ch] = p2_sx_2[ch];
				p2_sx_2[ch] = ret1; // inputs[ch][i];
				p2_sy_1[ch] = p2_sy_2[ch];
				p2_sy_2[ch] = ret2;
				double ret3 = c[ch][2] * (ret2 + p3_sy_1[ch]) - p3_sx_1[ch];
				p3_sx_1[ch] = p3_sx_2[ch];
				p3_sx_2[ch] = ret2; // inputs[ch][i];
				p3_sy_1[ch] = p3_sy_2[ch];
				p3_sy_2[ch] = ret3;
				double ret4 = c[ch][3] * (ret3 + p4_sy_1[ch]) - p4_sx_1[ch];
				p4_sx_1[ch] = p4_sx_2[ch];
				p4_sx_2[ch] = ret3; // inputs[ch][i];
				p4_sy_1[ch] = p4_sy_2[ch];
				p4_sy_2[ch] = ret4;

				if (ch == 1) outputs[ch][i] = p4_sy_1[ch];
				else outputs[ch][i] = p4_sy_2[ch];
			}
		}

		if (numChannels > 1) {
			for (int i = 0; i < sampleFrames; i++) {

				int ch = 1;
				double env = 0.637 * sqrt(outputs[0][i] * outputs[0][i] + outputs[1][i] * outputs[1][i]);
				double lp =
					(0.000042443309351420545) * env +
					(0.00008488661870284109) * sx_2[ch] +
					(0.000042443309351420545) * sx_1[ch] -
					(-1.9814857645620922) * sy_2[ch] -
					(0.9816555377994975) * sy_1[ch];
				sx_1[ch] = sx_2[ch];
				sx_2[ch] = env;
				sy_1[ch] = sy_2[ch];
				sy_2[ch] = lp;

				pp = (pp * coeff) + (icoef * env);

				// coef from 20ms ~ 2000ms
				rms_s = (rms_s * coeff) + (icoef * outputs[0][i] * outputs[0][i]);
				rms_c = (rms_c * coeff) + (icoef * outputs[1][i] * outputs[1][i]);
				double rhms = sqrt(rms_s + rms_c);

				double abs_in = std::abs(qq[0].front());
				rms = (rms * coeff) + (icoef * abs_in * abs_in);

				peak = (peak * coeff) + (icoef * abs_in);

				//outputs[1][i] = lp;  //Hilbert Low-passed
				//outputs[1][i] = env; //Hilbert as-is
				outputs[1][i] = rhms;
				//outputs[0][i] = qq[0].front(); // Original 
				//outputs[0][i] = sqrt(rms);
				outputs[0][i] = peak;
				//outputs[1][i] = pp;
				qq[0].pop();
				qq[1].pop();
			}
		}

		return;
	}
	
	/*
	Phase reference path c coefficients :
	0.47944111608296202665, 0.87624358989504858020, 0.97660296916871658368, 0.99749940412203375040,
	+90 deg path c coefficients :
	0.16177741706363166219, 0.73306690130335572242, 0.94536301966806279840, 0.99060051416704042460,
	*/


	compressor.setPower(bBypass);
	compressor.setThreshold(fThreshold * (0.0 - (-60.0)) + (-60.0));
	compressor.setRatio(fRatio * (20.0 - (1.0)) + (1.0));
	compressor.setKnee(fKnee * 20.0);
	compressor.setAttack(0.1 * exp(Attack_LOG_MAX * fAttack));
	compressor.setRelease(5.0 * exp(Release_LOG_MAX * fRelease));
	compressor.setLookAhead(fLookAhead * (512.0 - 0.0) + 0.0);
	compressor.setLookAheadEnable(bLookAheadEnable);
	compressor.setMakeup(fMakeup * (12.0 - (-12.0)) + (-12.0));
	compressor.setMix(fMix);

	//Update input peak metering
	inLevelFollower.updatePeak(inputs, numChannels, sampleFrames);


	/*
	// check setupProcessing
	for (int ch = 0; ch < numChannels; ch++) {
		for (int i = 0; i < sampleFrames; i++) {
			// OS - upsample
			// double in[2] = { 0.0, };
			// memmove(up[ch].buff + 1, up[ch].buff, sizeof(double) * (maxTap)-1);
			up[ch].acc();
			up[ch].buff[up[ch].now] = (double)inputs[ch][i];
			double in0 = 0.0;
			for (int i = 0; i < fir_size; i++) {
				//in0 += 2.0 * impulse_20240402[i] * up[ch].buff[up[ch].get_nth(i)];
				in0 += 2.0 * up[ch].coef[i] * up[ch].buff[up[ch].get_nth(i)];
			}
			//memmove(up[ch].buff + 1, up[ch].buff, sizeof(double) * (maxTap)-1);
			up[ch].acc();
			up[ch].buff[up[ch].now] = 0.0;
			double in1 = 0.0;
			for (int i = 0; i < fir_size; i++) {
				//in1 += 2.0 * impulse_20240402[i] * up[ch].buff[up[ch].get_nth(i)];
				in1 += 2.0 * up[ch].coef[i] * up[ch].buff[up[ch].get_nth(i)];
			}

			tt[ch][(2*i) +0] = in0;
			tt[ch][(2*i) +1] = in1;
		}
	}
	*/
	

	// Do compressor processing	
	compressor.process(\
		inputs, \
		outputs,\
		numChannels,\
		getSampleRate,\
		sampleFrames); // 2 * sampleFrames
		

	//memmove(ddd[0], ttt[0], sizeof(double) * 2 * sampleFrames);
	//memmove(ddd[1], ttt[1], sizeof(double) * 2 * sampleFrames);

	/*
	for (int ch = 0; ch < numChannels; ch++) {
		for (int i = 0; i < sampleFrames; i++) {
			// OS - downsample
			//memmove(dn[ch].buff + 1, dn[ch].buff, sizeof(double) * (maxTap)-1);
			dn[ch].acc();
			dn[ch].buff[dn[ch].now] = dd[ch][(2 * i) + 0];

			double out = 0.0;
			for (int i = 0; i < fir_size; i++) {
				//out += impulse_20240402[i] * dn[ch].buff[dn[ch].get_nth(i)];
				out += dn[ch].coef[i] * dn[ch].buff[dn[ch].get_nth(i)];
			}

			//memmove(dn[ch].buff + 1, dn[ch].buff, sizeof(double) * (maxTap)-1);
			dn[ch].acc();
			dn[ch].buff[dn[ch].now] = dd[ch][(2 * i) + 1];

			outputs[ch][i] = (SampleType)out;
		}
	}
	*/

	// Update gain reduction metering
	gainReduction = compressor.getMaxGainReduction();
	
	// Update output peak metering
	outLevelFollower.updatePeak(outputs, numChannels, sampleFrames);


	for (int ch = 0; ch < numChannels; ch++) {
		currentInput[ch] = inLevelFollower.getPeak(ch);
		currentOutput[ch] = outLevelFollower.getPeak(ch);
	}
	
	if (currentExchangeBlock.blockID == Vst::InvalidDataExchangeBlockID)
		acquireNewExchangeBlock();

	if (auto block = toDataBlock(currentExchangeBlock))
	{
		auto numSamples = static_cast<uint32> (sampleFrames);
		while (numSamples > 0)
		{
			uint32 numSamplesFreeInBlock = block->sampleRate - block->numSamples;
			uint32 numSamplesToCopy = std::min<uint32>(numSamplesFreeInBlock, numSamples);

			/*
			for (auto channel = 0; channel < input.numChannels; ++channel)
			{
				auto blockChannelData = &block->samples[0] + block->numSamples;
				auto inputChannel = input.channelBuffers32[channel] + (data.numSamples - numSamples);
				memcpy(blockChannelData, inputChannel, numSamplesToCopy * sizeof(float));
			}
			*/
			block->numSamples += numSamplesToCopy;

			if (block->numSamples * 30 >= block->sampleRate) // 1sec
			{
				block->inL = currentInput[0];
				block->inR = currentInput[1];
				block->outL = currentOutput[0];
				block->outR = currentOutput[1];
				block->gR = gainReduction;

				dataExchange->sendCurrentBlock();
				acquireNewExchangeBlock();
				block = toDataBlock(currentExchangeBlock);
				if (block == nullptr)
					break;
			}

			numSamples -= numSamplesToCopy;
		}
	}

	/*
	cntt += sampleFrames;
	if (cntt * 30 > sampleFrames) {
		for (int ch = 0; ch < numChannels; ch++) {
			currentInput[ch] = Decibels::gainToDecibels(inLevelFollower.getPeak(ch));
			auto msg = owned(allocateMessage());
			msg->setMessageID("In");
			msg->getAttributes()->setFloat(std::to_string(ch).c_str(), currentInput[ch]);
			sendMessage(msg);
		}
		for (int ch = 0; ch < numChannels; ch++) {
			currentOutput[ch] = Decibels::gainToDecibels(outLevelFollower.getPeak(ch));
			auto msg = owned(allocateMessage());
			msg->setMessageID("Out");
			msg->getAttributes()->setFloat(std::to_string(ch).c_str(), currentOutput[ch]);
			sendMessage(msg);
		}
		{
			auto msg = owned(allocateMessage());
			msg->setMessageID("GR");
			msg->getAttributes()->setFloat(std::to_string(0).c_str(), gainReduction);
			sendMessage(msg);
		}
		cntt = 0;
	}
	*/


	return;
}

//------------------------------------------------------------------------
} // namespace yg331