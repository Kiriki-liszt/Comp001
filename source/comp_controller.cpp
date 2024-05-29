//------------------------------------------------------------------------
// Copyright(c) 2024 yg331.
//------------------------------------------------------------------------

#include "comp_controller.h"
//#include "comp_cids.h"
#include "vstgui/plugin-bindings/vst3editor.h"
#include "pluginterfaces/base/ustring.h"
#include "base/source/fstreamer.h"
#include "base/source/fdebug.h"

#include "vstgui/vstgui.h"
#include "vstgui/vstgui_uidescription.h"
#include "vstgui/uidescription/detail/uiviewcreatorattributes.h"
#include "vstgui/lib/controls/cvumeter.h"
#include "vstgui/lib/controls/cparamdisplay.h"

static const std::string kAttrBackColor = "back-color";
static const std::string kAttrLineColor = "line-color";
static const std::string kAttrVuOnColor = "vu-on-color";
static const std::string kAttrVuOffColor = "vu-off-color";
static const std::string kAttrPDclick = "click-behave";
static const std::string kMin = "Min";
static const std::string kMax = "Max";
using namespace Steinberg;
namespace VSTGUI {
	class MyControlFactory : public ViewCreatorAdapter
	{
	public:
		//register this class with the view factory
		MyControlFactory() { UIViewFactory::registerViewCreator(*this); }

		//return an unique name here
		IdStringPtr getViewName() const override { return "My Vu Meter"; }

		//return the name here from where your custom view inherites.
		//	Your view automatically supports the attributes from it.
		IdStringPtr getBaseViewName() const override { return UIViewCreator::kCControl; }

		//create your view here.
		//	Note you don't need to apply attributes here as
		//	the apply method will be called with this new view
		CView* create(const UIAttributes & attributes, const IUIDescription * description) const override
		{
			CRect size(CPoint(45, 45), CPoint(400, 150));
			return new MyVuMeter(size, 2);
		}
		bool apply(
			CView* view, 
			const UIAttributes& attributes,
			const IUIDescription* description) const
		{
			auto* vuMeter = dynamic_cast<MyVuMeter*> (view);

			if (!vuMeter)
				return false;

			const auto* attr = attributes.getAttributeValue(UIViewCreator::kAttrOrientation);
			if (attr)
				vuMeter->setStyle(*attr == UIViewCreator::strVertical ? MyVuMeter::kVertical : MyVuMeter::kHorizontal);

			CColor color;
			if (UIViewCreator::stringToColor(attributes.getAttributeValue(kAttrVuOnColor), color, description))
				vuMeter->setVuOnColor(color);
			if (UIViewCreator::stringToColor(attributes.getAttributeValue(kAttrVuOffColor), color, description))
				vuMeter->setVuOffColor(color);

			return true;
		}

		bool getAttributeNames(StringList& attributeNames) const
		{
			attributeNames.emplace_back(UIViewCreator::kAttrOrientation);
			attributeNames.emplace_back(kAttrVuOnColor);
			attributeNames.emplace_back(kAttrVuOffColor);
			return true;
		}

		AttrType getAttributeType(const std::string& attributeName) const
		{
			if (attributeName == UIViewCreator::kAttrOrientation)
				return kListType;
			if (attributeName == kAttrVuOnColor)
				return kColorType;
			if (attributeName == kAttrVuOffColor)
				return kColorType;
			return kUnknownType;
		}

		//------------------------------------------------------------------------
		bool getAttributeValue(
			CView* view, 
			const string& attributeName,
			string& stringValue, 
			const IUIDescription* desc) const
		{
			auto* vuMeter = dynamic_cast<MyVuMeter*> (view);

			if (!vuMeter)
				return false;

			if (attributeName == UIViewCreator::kAttrOrientation)
			{
				if (vuMeter->getStyle() & MyVuMeter::kVertical)
					stringValue = UIViewCreator::strVertical;
				else
					stringValue = UIViewCreator::strHorizontal;
				return true;
			}
			else if (attributeName == kAttrVuOnColor)
			{
				UIViewCreator::colorToString(vuMeter->getVuOnColor(), stringValue, desc);
				return true;
			}
			else if (attributeName == kAttrVuOffColor)
			{
				UIViewCreator::colorToString(vuMeter->getVuOffColor(), stringValue, desc);
				return true;
			}
			return false;
		}

		//------------------------------------------------------------------------
		bool getPossibleListValues(
			const string& attributeName,
			ConstStringPtrList& values) const
		{
			if (attributeName == UIViewCreator::kAttrOrientation)
			{
				return UIViewCreator::getStandardAttributeListValues(UIViewCreator::kAttrOrientation, values);
			}
			return false;
		}

	};





	class MyCurveControlFactory : public ViewCreatorAdapter
	{
	public:
		MyCurveControlFactory() { UIViewFactory::registerViewCreator(*this); }
		IdStringPtr getViewName() const override { return "My Curve Viewer"; }
		IdStringPtr getBaseViewName() const override { return UIViewCreator::kCControl; }
		CView* create(const UIAttributes& attributes, const IUIDescription* description) const override
		{
			CRect size(CPoint(45, 45), CPoint(400, 150));
			return new CurveView(size);
		}
		bool apply(
			CView* view,
			const UIAttributes& attributes,
			const IUIDescription* description) const
		{
			auto* curveView = dynamic_cast<CurveView*> (view);

			if (!curveView)
				return false;

			CColor color;
			if (UIViewCreator::stringToColor(attributes.getAttributeValue(kAttrBackColor), color, description))
				curveView->setBackColor(color);
			if (UIViewCreator::stringToColor(attributes.getAttributeValue(kAttrLineColor), color, description))
				curveView->setLineColor(color);

			return true;
		}

		bool getAttributeNames(StringList& attributeNames) const
		{
			attributeNames.emplace_back(kAttrBackColor);
			attributeNames.emplace_back(kAttrLineColor);
			return true;
		}

		AttrType getAttributeType(const std::string& attributeName) const
		{
			if (attributeName == kAttrBackColor)
				return kColorType;
			if (attributeName == kAttrLineColor)
				return kColorType;
			return kUnknownType;
		}

		//------------------------------------------------------------------------
		bool getAttributeValue(
			CView* view,
			const string& attributeName,
			string& stringValue,
			const IUIDescription* desc) const
		{
			auto* curveView = dynamic_cast<CurveView*> (view);

			if (!curveView)
				return false;

			if (attributeName == kAttrVuOnColor)
			{
				UIViewCreator::colorToString(curveView->getBackColor(), stringValue, desc);
				return true;
			}
			else if (attributeName == kAttrVuOffColor)
			{
				UIViewCreator::colorToString(curveView->getLineColor(), stringValue, desc);
				return true;
			}
			return false;
		}
	};


	class MyPDFactory : public ViewCreatorAdapter
	{
	public:
		MyPDFactory() { UIViewFactory::registerViewCreator(*this); }
		IdStringPtr getViewName() const override { return "My PD"; }
		IdStringPtr getBaseViewName() const override { return UIViewCreator::kCParamDisplay; }
		CView* create(const UIAttributes& attributes, const IUIDescription* description) const override
		{
			CRect ss(0, 0, 100, 20);
			return new MyPD(ss, nullptr, 0);
		}
		bool apply(
			CView* view,
			const UIAttributes& attributes,
			const IUIDescription* description) const
		{
			auto* vv = dynamic_cast<MyPD*> (view);

			if (!vv)
				return false;

			const auto* attr = attributes.getAttributeValue(kAttrPDclick);
			if (attr)
				vv->setStyle_(*attr == kMin ? MyPD::kMin : MyPD::kMax);

			return true;
		}

		bool getAttributeNames(StringList& attributeNames) const
		{
			attributeNames.emplace_back(kAttrPDclick);
			return true;
		}

		AttrType getAttributeType(const std::string& attributeName) const
		{
			if (attributeName == kAttrPDclick)
				return kListType;
			return kUnknownType;
		}

		//------------------------------------------------------------------------
		bool getAttributeValue(
			CView* view,
			const string& attributeName,
			string& stringValue,
			const IUIDescription* desc) const
		{
			auto* vv = dynamic_cast<MyPD*> (view);

			if (!vv)
				return false;

			if (attributeName == kAttrPDclick)
			{
				if (vv->getStyle_() & MyPD::kMin)
					stringValue = kMin;
				else
					stringValue = kMax;
				return true;
			}

			return false;
		}

		//------------------------------------------------------------------------
		bool getPossibleListValues(
			const string& attributeName,
			ConstStringPtrList& values) const
		{
			if (attributeName == kAttrPDclick)
			{
				values.emplace_back(&kMin);
				values.emplace_back(&kMax);
				return true;
			}
			return false;
		}
	};

	class MeterViewContainerFactory : public ViewCreatorAdapter
	{
	public:
		MeterViewContainerFactory() { UIViewFactory::registerViewCreator(*this); }
		IdStringPtr getViewName() const override { return "MeterViewContainer"; }
		IdStringPtr getBaseViewName() const override { return UIViewCreator::kCViewContainer; }
		CView* create(const UIAttributes& attributes, const IUIDescription* description) const override
		{
			CRect ss(0, 0, 100, 20);
			return new MeterViewContainer(ss);
		}
		
	};

	//create a static instance so that it registers itself with the view factory
	MyControlFactory __gMyControlFactory;
	MyCurveControlFactory __gMyCurveControlFactory;
	MyPDFactory __gMyPDFactory;
	MeterViewContainerFactory __gMeterViewContainerFactory;
} // namespace VSTGUI

namespace yg331 {
//------------------------------------------------------------------------
// LogRangeParameter Declaration
//------------------------------------------------------------------------
	class LogRangeParameter : public Vst::RangeParameter
	{
	public:
		using RangeParameter::RangeParameter;
		Vst::ParamValue toPlain(Vst::ParamValue _valueNormalized) const SMTG_OVERRIDE;
		Vst::ParamValue toNormalized(Vst::ParamValue plainValue) const SMTG_OVERRIDE;
		void toString(Vst::ParamValue _valueNormalized, Vst::String128 string) const SMTG_OVERRIDE;
	};
	//------------------------------------------------------------------------
	// LogRangeParameter Implementation
	//------------------------------------------------------------------------
	Vst::ParamValue LogRangeParameter::toPlain(Vst::ParamValue _valueNormalized) const
	{
		double FREQ_LOG_MAX = log(getMax() / getMin());
		double tmp = getMin() * exp(FREQ_LOG_MAX * _valueNormalized);
		double freq = (std::max)((std::min)(tmp, getMax()), getMin());
		return freq;
		//return _valueNormalized * (getMax() - getMin()) + getMin();
	}

	//------------------------------------------------------------------------
	Vst::ParamValue LogRangeParameter::toNormalized(Vst::ParamValue plainValue) const
	{
		SMTG_ASSERT(getMax() - getMin() != 0);
		double FREQ_LOG_MAX = log(getMax() / getMin());
		return log(plainValue / getMin()) / FREQ_LOG_MAX;
		//return (plainValue - getMin()) / (getMax() - getMin());
	}

	void LogRangeParameter::toString(Vst::ParamValue _valueNormalized, Vst::String128 string) const
	{
		{
			//Parameter::toString(toPlain(_valueNormalized), string);
			UString wrapper(string, str16BufferSize(Vst::String128));
			{
				if (!wrapper.printFloat(toPlain(_valueNormalized), precision))
					string[0] = 0;
				wrapper.append(STR16(" "));
				wrapper.append(getInfo().units);
			}
		}
	}

	//------------------------------------------------------------------------
	// LinRangeParameter Declaration
	//------------------------------------------------------------------------
	class LinRangeParameter : public Vst::RangeParameter
	{
	public:
		using RangeParameter::RangeParameter;
		void toString(Vst::ParamValue _valueNormalized, Vst::String128 string) const SMTG_OVERRIDE;
	};
	//------------------------------------------------------------------------
	// LinRangeParameter Implementation
	//------------------------------------------------------------------------
	void LinRangeParameter::toString(Vst::ParamValue _valueNormalized, Vst::String128 string) const
	{
		{
			//Parameter::toString(toPlain(_valueNormalized), string);
			UString wrapper(string, str16BufferSize(Vst::String128));
			{
				if (!wrapper.printFloat(toPlain(_valueNormalized), precision))
					string[0] = 0;
				wrapper.append(STR16(" "));
				wrapper.append(getInfo().units);
			}
		}
	}
//------------------------------------------------------------------------
// comp_Controller Implementation
//------------------------------------------------------------------------
tresult PLUGIN_API comp_Controller::initialize (FUnknown* context)
{
	// Here the Plug-in will be instantiated

	//---do not forget to call parent ------
	tresult result = EditControllerEx1::initialize (context);
	if (result != kResultOk)
	{
		return result;
	}

	// Here you could register some parameters

	int32 stepCount;
	int32 flags;
	int32 tag;
	Vst::ParamValue defaultVal;
	Vst::ParamValue defaultPlain;


	tag = kParamBypass;
	stepCount = 1;
	defaultVal = Init_Bypass ? 1 : 0;
	flags = Vst::ParameterInfo::kCanAutomate | Vst::ParameterInfo::kIsBypass;
	parameters.addParameter(STR16("Bypass"), nullptr, stepCount, defaultVal, flags, tag);

	if (zoomFactors.empty())
	{
		Vst::ParamValue zoom_coef = 1.0;
		zoomFactors.push_back(ZoomFactor(STR("50%"), zoom_coef * 0.5));  // 0/6
		zoomFactors.push_back(ZoomFactor(STR("75%"), zoom_coef * 0.75)); // 1/6
		zoomFactors.push_back(ZoomFactor(STR("100%"), zoom_coef * 1.0));  // 2/6
		zoomFactors.push_back(ZoomFactor(STR("125%"), zoom_coef * 1.25)); // 3/6
		zoomFactors.push_back(ZoomFactor(STR("150%"), zoom_coef * 1.5));  // 4/6
		zoomFactors.push_back(ZoomFactor(STR("175%"), zoom_coef * 1.75)); // 5/6
		zoomFactors.push_back(ZoomFactor(STR("200%"), zoom_coef * 2.0));  // 6/6
	}

	auto zoomParameter = new Vst::StringListParameter(STR("Zoom"), kParamZoom);
	for (auto it = zoomFactors.begin(), end = zoomFactors.end(); it != end; ++it)
	{
		zoomParameter->appendString(it->title);
	}
	zoomParameter->setNormalized(zoomParameter->toNormalized(Init_Zoom)); // toNorm(2) == 100%
	zoomParameter->addDependent(this);
	parameters.addParameter(zoomParameter);


	flags = Vst::ParameterInfo::kCanAutomate;
	defaultPlain = 0.0;
	stepCount = 0;

	tag = kParamThreshold;
	auto* ParamThreshold = new LinRangeParameter(STR16("Threshold"), tag, STR16("dB"), -60, 0, 0, stepCount, flags);
	ParamThreshold->setPrecision(2);
	parameters.addParameter(ParamThreshold);

	tag = kParamRatio;
	auto* ParamRatio = new Vst::RangeParameter(STR16("Ratio"), tag, STR16(""), 1, 20, 4, stepCount, flags);
	ParamRatio->setPrecision(2);
	parameters.addParameter(ParamRatio);

	tag = kParamKnee;
	auto* ParamKnee = new LinRangeParameter(STR16("Knee"), tag, STR16("dB"), 0, 20, 5, stepCount, flags);
	ParamKnee->setPrecision(2);
	parameters.addParameter(ParamKnee);

	tag = kParamAttack;
	auto* ParamAttack = new LogRangeParameter(STR16("Attack"), tag, STR16("ms"), 0.1, 50, 1, stepCount, flags);
	ParamAttack->setPrecision(2);
	parameters.addParameter(ParamAttack);

	tag = kParamRelease;
	auto* ParamRelease = new LogRangeParameter(STR16("Release"), tag, STR16("ms"), 5.0, 5000.0, 50.0, stepCount, flags);
	ParamRelease->setPrecision(2);
	parameters.addParameter(ParamRelease);
	
	tag = kParamLookAheadEnable;
	parameters.addParameter(STR16("LookAhead_Enable"), nullptr, 1 /*stepCount*/, Init_Bypass ? 1 : 0/*defaultVal*/, flags, tag);

	tag = kParamLookAhead;
	auto* ParamLookAhead = new Vst::RangeParameter(STR16("LookAhead"), tag, STR16("smpls"), 0, 512, 0, stepCount, flags);
	ParamLookAhead->setPrecision(2);
	parameters.addParameter(ParamLookAhead);

	tag = kParamMakeup;
	auto* ParamMakeup = new LinRangeParameter(STR16("Makeup"), tag, STR16("dB"), -12, 12, 0, stepCount, flags);
	ParamMakeup->setPrecision(2);
	parameters.addParameter(ParamMakeup);

	tag = kParamMix;
	auto* ParamMix = new LinRangeParameter(STR16("Mix"), tag, STR16("%"), 0, 100, 100, stepCount, flags);
	ParamMix->setPrecision(2);
	parameters.addParameter(ParamMix);
	
	return result;
}

//------------------------------------------------------------------------
tresult PLUGIN_API comp_Controller::terminate ()
{
	// Here the Plug-in will be de-instantiated, last possibility to remove some memory!

	//---do not forget to call parent ------
	return EditControllerEx1::terminate ();
}

//------------------------------------------------------------------------
tresult PLUGIN_API comp_Controller::setComponentState (IBStream* state)
{
	// Here you get the state of the component (Processor part)
	if (!state)
		return kResultFalse;

	IBStreamer streamer(state, kLittleEndian);

	int32           savedBypass          = 0;
	Vst::ParamValue savedZoom            = 0.0;

	Vst::ParamValue savedThreshold       = 0.0;
	Vst::ParamValue savedRatio           = 0.0;
	Vst::ParamValue savedKnee            = 0.0;
	Vst::ParamValue savedAttack          = 0.0;
	Vst::ParamValue savedRelease         = 0.0;
	Vst::ParamValue savedLookAhead       = 0.0;
	int32           savedLookAheadEnable = 0;

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

	setParamNormalized(kParamBypass,          savedBypass ? 1 : 0);
	setParamNormalized(kParamZoom,            savedZoom);
	setParamNormalized(kParamThreshold,       savedThreshold);
	setParamNormalized(kParamRatio,           savedRatio);
	setParamNormalized(kParamKnee,            savedKnee);
	setParamNormalized(kParamAttack,          savedAttack);
	setParamNormalized(kParamRelease,         savedRelease);
	setParamNormalized(kParamLookAhead,       savedLookAhead);
	setParamNormalized(kParamLookAheadEnable, savedLookAheadEnable ? 1 : 0);
	setParamNormalized(kParamMakeup,          savedMakeup);
	setParamNormalized(kParamMix,             savedMix);

	return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API comp_Controller::setState (IBStream* state)
{
	// Here you get the state of the controller

	return kResultTrue;
}

//------------------------------------------------------------------------
tresult PLUGIN_API comp_Controller::getState (IBStream* state)
{
	// Here you are asked to deliver the state of the controller (if needed)
	// Note: the real state of your plug-in is saved in the processor

	return kResultTrue;
}

//------------------------------------------------------------------------
IPlugView* PLUGIN_API comp_Controller::createView (FIDString name)
{
	// Here the Host wants to open your editor (if you have one)
	if (FIDStringsEqual (name, Vst::ViewType::kEditor))
	{
		// create your editor here and return a IPlugView ptr of it
		auto* view = new VSTGUI::VST3Editor (this, "view", "comp_editor.uidesc");
		view->setZoomFactor(1.0);
		setKnobMode(Steinberg::Vst::KnobModes::kLinearMode);
		eded = view;
		return view;
	}
	return nullptr;
}

void PLUGIN_API comp_Controller::update(FUnknown* changedUnknown, int32 message)
{
	EditControllerEx1::update(changedUnknown, message);

	// GUI Resizing
	// check 'zoomtest' code at
	// https://github.com/steinbergmedia/vstgui/tree/vstgui4_10/vstgui/tests/uidescription%20vst3/source

	Vst::Parameter* param = FCast<Vst::Parameter>(changedUnknown);
	if (!param)
		return;

	if (param->getInfo().id == kParamZoom)
	{
		size_t index = static_cast<size_t> (param->toPlain(param->getNormalized()));

		if (index >= zoomFactors.size())
			return;

		if (eded)
			eded->setZoomFactor(zoomFactors[index].factor);

		/*
		for (EditorVector::const_iterator it = editors.begin(), end = editors.end(); it != end; ++it)
		{
			VSTGUI::VST3Editor* editor = dynamic_cast<VSTGUI::VST3Editor*>(*it);
			if (editor)
				editor->setZoomFactor(zoomFactors[index].factor);
		}
		*/
	}
}
/*
//------------------------------------------------------------------------
void comp_Controller::editorAttached(Steinberg::Vst::EditorView* editor)
{
	editors.push_back(editor);
}

//------------------------------------------------------------------------
void comp_Controller::editorRemoved(Steinberg::Vst::EditorView* editor)
{
	editors.erase(std::find(editors.begin(), editors.end(), editor));
}
*/
//------------------------------------------------------------------------
tresult PLUGIN_API comp_Controller::setParamNormalized (Vst::ParamID tag, Vst::ParamValue value)
{
	// called by host to update your parameters
	tresult result = EditControllerEx1::setParamNormalized (tag, value);
	return result;
}

//------------------------------------------------------------------------
tresult PLUGIN_API comp_Controller::getParamStringByValue (Vst::ParamID tag, Vst::ParamValue valueNormalized, Vst::String128 string)
{
	// called by host to get a string for given normalized value of a specific parameter
	// (without having to set the value!)
	return EditControllerEx1::getParamStringByValue (tag, valueNormalized, string);
}

//------------------------------------------------------------------------
tresult PLUGIN_API comp_Controller::getParamValueByString (Vst::ParamID tag, Vst::TChar* string, Vst::ParamValue& valueNormalized)
{
	// called by host to get a normalized value from a string representation of a specific parameter
	// (without having to set the value!)
	return EditControllerEx1::getParamValueByString (tag, string, valueNormalized);
}

//------------------------------------------------------------------------
// DataExchangeController Implementation
//------------------------------------------------------------------------
tresult PLUGIN_API comp_Controller::notify(Vst::IMessage* message)
{
	if (dataExchange.onMessage(message))
		return kResultTrue;


	if (!message)
		return kInvalidArgument;
	
	/*
	if (FIDStringsEqual(message->getMessageID(), "In"))
	{
		double In = 0.0;
		if (message->getAttributes()->getFloat("0", In) == kResultOk)
		{
			InMeter[0] = In;
			if (!uiMessageControllers.empty()) {
				for (auto iter = uiMessageControllers.begin(); iter != uiMessageControllers.end(); iter++) {
					//	(*iter)->setValueDD(InMeter[0]);
					(*iter)->setVuMeterValue(InMeter[0], kInL);
				}
			}
		}
		if (message->getAttributes()->getFloat("1", In) == kResultOk)
		{
			InMeter[1] = In;
			if (!uiMessageControllers.empty()) {
				for (auto iter = uiMessageControllers.begin(); iter != uiMessageControllers.end(); iter++) {
					// (*iter)->setValueDD(InMeter[1]);
					(*iter)->setVuMeterValue(InMeter[1], kInR);
				}
			}
		}
		message->release();
	}

	if (FIDStringsEqual(message->getMessageID(), "Out"))
	{
		double In = 0.0;
		if (message->getAttributes()->getFloat("0", In) == kResultOk)
		{
			OutMeter[0] = In;
			if (!uiMessageControllers.empty()) {
				for (auto iter = uiMessageControllers.begin(); iter != uiMessageControllers.end(); iter++) {
					//	(*iter)->setValueDD(InMeter[0]);
					(*iter)->setVuMeterValue(OutMeter[0], kOutL);
				}
			}
		}
		if (message->getAttributes()->getFloat("1", In) == kResultOk)
		{
			OutMeter[1] = In;
			if (!uiMessageControllers.empty()) {
				for (auto iter = uiMessageControllers.begin(); iter != uiMessageControllers.end(); iter++) {
					// (*iter)->setValueDD(InMeter[1]);
					(*iter)->setVuMeterValue(OutMeter[1], kOutR);
				}
			}
		}
		message->release();
	}

	if (FIDStringsEqual(message->getMessageID(), "GR"))
	{
		double In = 0.0;
		if (message->getAttributes()->getFloat("0", In) == kResultOk)
		{
			GRMeter = In;
			if (!uiMessageControllers.empty()) {
				for (auto iter = uiMessageControllers.begin(); iter != uiMessageControllers.end(); iter++) {
					//	(*iter)->setValueDD(InMeter[0]);
					(*iter)->setVuMeterValue(GRMeter, kGR);
				}
			}
		}
		message->release();
	}
	*/
	
	return EditControllerEx1::notify(message);
}


//------------------------------------------------------------------------
void PLUGIN_API comp_Controller::queueOpened(Vst::DataExchangeUserContextID userContextID,
	uint32 blockSize,
	TBool& dispatchOnBackgroundThread)
{
	//FDebugPrint("Data Exchange Queue opened.\n");
}

//------------------------------------------------------------------------
void PLUGIN_API comp_Controller::queueClosed(Vst::DataExchangeUserContextID userContextID)
{
	//FDebugPrint("Data Exchange Queue closed.\n");
}

//------------------------------------------------------------------------
void PLUGIN_API comp_Controller::onDataExchangeBlocksReceived(
	Vst::DataExchangeUserContextID userContextID, 
	uint32 numBlocks, 
	Vst::DataExchangeBlock* blocks,
	TBool onBackgroundThread
)
{
	for (auto index = 0u; index < numBlocks; ++index)
	{
		auto dataBlock = toDataBlock(blocks[index]);
		InMeter[0] = dataBlock->inL;
		InMeter[1] = dataBlock->inR;
		OutMeter[0] = dataBlock->outL;
		OutMeter[1] = dataBlock->outR;
		GRMeter = dataBlock->gR;
		if (!vuMeterControllers.empty()) {
			for (auto iter = vuMeterControllers.begin(); iter != vuMeterControllers.end(); iter++) {
				(*iter)->setMeterValue(InMeter[0]*0.5 + InMeter[1]*0.5, kIn);
				(*iter)->setMeterValue(OutMeter[0]*0.5 + OutMeter[1]*0.5, kOut);
				(*iter)->setMeterValue(GRMeter, kGR);
				(*iter)->setVuMeterValue(
					InMeter[0], InMeter[1],
					OutMeter[0], OutMeter[1],
					GRMeter
				);
				/*
				(*iter)->setVuMeterValue(InMeter[0], kInL);
				(*iter)->setVuMeterValue(InMeter[1], kInR);
				(*iter)->setVuMeterValue(OutMeter[0], kOutL);
				(*iter)->setVuMeterValue(OutMeter[1], kOutR);
				(*iter)->setVuMeterValue(GRMeter, kGR);
				*/
			}
		}
		/*
		FDebugPrint(
			"Received Data Block: SampleRate: %d, SampleSize: %d, NumChannels: %d, NumSamples: %d\n",
			dataBlock->sampleRate, 
			static_cast<uint32_t> (dataBlock->sampleSize),
			static_cast<uint32_t> (dataBlock->numChannels),
			static_cast<uint32_t> (dataBlock->numSamples));
		*/
		//FDebugPrint(\
			"Received Data Block: %f %f %f %f %f\n",\
			dataBlock->inL,\
			dataBlock->inR,\
			dataBlock->outL,\
			dataBlock->outR,\
			dataBlock->gR);
	}
}

//------------------------------------------------------------------------
} // namespace yg331
