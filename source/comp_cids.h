//------------------------------------------------------------------------
// Copyright(c) 2024 yg331.
//------------------------------------------------------------------------

#pragma once

#include "pluginterfaces/base/funknown.h"
#include "pluginterfaces/vst/vsttypes.h"

enum {
	kParamBypass = 0,
	kParamZoom,

	kParamThreshold,
	kParamRatio,
	kParamKnee,

	kParamAttack,
	kParamRelease,
	kParamLookAhead,
	kParamLookAheadEnable,

	kParamMakeup,
	kParamMix
};

enum Tags
{
	kIn = 100,
	kInL = 101,
	kInR = 102,
	kOut = 103,
	kOutL = 104,
	kOutR = 105,
	kGR = 106,
	
};

const bool
Init_Bypass = false,
Init_LookAheadEnable = false;

const Steinberg::Vst::ParamValue
Init_Zoom = 2.0;

namespace yg331 {
//------------------------------------------------------------------------
static const Steinberg::FUID kcomp_ProcessorUID (0x973A8359, 0x73F55212, 0xBA5BA547, 0x6DAC162B);
static const Steinberg::FUID kcomp_ControllerUID (0xD2C4FDC6, 0x542D5371, 0xBEC41EB9, 0x38816A52);

#define comp_VST3Category "Fx"

//------------------------------------------------------------------------
} // namespace yg331
