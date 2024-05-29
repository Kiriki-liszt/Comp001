//------------------------------------------------------------------------
// Copyright(c) 2024 yg331.
//------------------------------------------------------------------------

#pragma once
#include "comp_cids.h"
#include "comp_dataexchange.h"
#include "public.sdk/source/vst/vsteditcontroller.h"
#include "vstgui/plugin-bindings/vst3editor.h"
#include "vstgui/uidescription/delegationcontroller.h"
#include "vstgui/lib/iviewlistener.h"
#include "vstgui/uidescription/icontroller.h"
#include "public.sdk/source/vst/utility/stringconvert.h"
#include "base/source/fstring.h"

using namespace Steinberg;

namespace VSTGUI {
	//------------------------------------------------------------------------
	//  Parameter Display for Min/Max indication
	//------------------------------------------------------------------------
	class MyPD : public CParamDisplay {
	private:
		enum StyleEnum
		{
			StyleMin = 0,
			StyleMax,
		};
	public:
		enum Style
		{
			kMin = 1 << StyleMin,
			kMax = 1 << StyleMax,
		};
		MyPD(const CRect& size, CBitmap* background, int32_t inStyle, int32_t style = kMin)
			: CParamDisplay(size, background, inStyle), _style(style)
		{};

		MyPD(const CParamDisplay& v, int32_t style = kMin)
			: CParamDisplay(v), _style(style)
		{};
		void setValue(float val) override {
			double vv = 0.0;
			if (_style == kMax) 
				vv = std::max(getValue(), val);
			else 
				vv = std::min(getValue(), val);
			CParamDisplay::setValue(vv);
		};
		void setStyle_(int32_t newStyle) { _style = newStyle; invalid(); }
		int32_t getStyle_() const { return _style; }
		int32_t     _style;
	};


	//------------------------------------------------------------------------
	//  VU meter view
	//------------------------------------------------------------------------
	class MyVuMeter : public CControl {
	private:
		enum StyleEnum
		{
			StyleHorizontal = 0,
			StyleVertical,
		};
	public:
		enum Style
		{
			kHorizontal = 1 << StyleHorizontal,
			kVertical = 1 << StyleVertical,
		};

		MyVuMeter(const CRect& size, int32_t style = kVertical) 
			: CControl(size, nullptr, 0)
			, style(style)
		{
			vuOnColor = kWhiteCColor;
			vuOffColor = kBlackCColor;

			rectOn(size.left, size.top, size.right, size.bottom);
			rectOff(size.left, size.top, size.right, size.bottom);

			setWantsIdle(true);
		}
		MyVuMeter(const MyVuMeter& vuMeter)
			: CControl(vuMeter)
			, style(vuMeter.style)
			, vuOnColor(vuMeter.vuOnColor)
			, vuOffColor(vuMeter.vuOffColor)
			, rectOn(vuMeter.rectOn)
			, rectOff(vuMeter.rectOff)
		{
			setWantsIdle(true);
		}

		void setStyle(int32_t newStyle) { style = newStyle; invalid(); }
		int32_t getStyle() const { return style; }

		virtual void setVuOnColor(CColor color) {
			if (vuOnColor != color) { vuOnColor = color; setDirty(true); }
		}
		CColor getVuOnColor() const { return vuOnColor; }

		virtual void setVuOffColor(CColor color) {
			if (vuOffColor != color) { vuOffColor = color; setDirty(true); }
		}
		CColor getVuOffColor() const { return vuOffColor; }

		// overrides
		void setDirty(bool state) override
		{
			CView::setDirty(state);
		};
		void draw(CDrawContext* _pContext) override {

			CRect _rectOn(rectOn);
			CRect _rectOff(rectOff);
			CPoint pointOn;
			CPoint pointOff;
			CDrawContext* pContext = _pContext;

			bounceValue();

			float newValue = getValueNormalized(); // normalize

			if (style & kHorizontal)
			{
				auto tmp = (CCoord)((int32_t)(newValue) * getViewSize().getWidth());
				pointOff(tmp, 0);

				_rectOff.left += tmp;
				_rectOn.right = tmp + rectOn.left;
			}
			else
			{
				auto tmp = (CCoord)((int32_t)(newValue * getViewSize().getHeight()));
				pointOn(0, tmp);

				//_rectOff.bottom = tmp + rectOff.top;
				//_rectOn.top += tmp;
				_rectOn.top = _rectOff.bottom - tmp;
			}

			pContext->setFillColor(vuOffColor);
			pContext->drawRect(rectOff, kDrawFilled);

			pContext->setFillColor(vuOnColor);
			pContext->drawRect(_rectOn, kDrawFilled);

			setDirty(false);
		};
		void setViewSize(const CRect& newSize, bool invalid = true) override
		{
			CControl::setViewSize(newSize, invalid);
			rectOn = getViewSize();
			rectOff = getViewSize();
		};
		bool sizeToFit() override {
			if (getDrawBackground())
			{
				CRect vs(getViewSize());
				vs.setWidth(getDrawBackground()->getWidth());
				vs.setHeight(getDrawBackground()->getHeight());
				setViewSize(vs);
				setMouseableArea(vs);
				return true;
			}
			return false;
		};
		/*
		void onIdle() override {
			if (getOldValue() != value)
				invalid();
		};
		*/

		CLASS_METHODS(MyVuMeter, CControl)

	protected:
		~MyVuMeter() noexcept override
		{
			//setOnBitmap(nullptr);
			//setOffBitmap(nullptr);
		};

		int32_t     style;

		CColor		vuOnColor;
		CColor		vuOffColor;

		CRect    rectOn;
		CRect    rectOff;
	};



	//------------------------------------------------------------------------
	//  Curve Display
	//------------------------------------------------------------------------
	class CurveView : public CControl {
	public:
		CurveView(
			const VSTGUI::CRect& size, 
			Steinberg::Vst::EditController* _editController = nullptr)
			: CControl(size, nullptr, 0)
		{
			BackColor = kWhiteCColor;
			LineColor = kBlackCColor;
			setWantsIdle(true);
		}
		CurveView(const CurveView& vuMeter)
			: CControl(vuMeter)
			, BackColor(vuMeter.BackColor)
			, LineColor(vuMeter.LineColor)
		{
			setWantsIdle(true);
		}

		// get/set Parameters
		virtual void setThreshold(double value) { if (Threshold != value) { Threshold = value; setDirty(true); } }
		double getThrershold() const { return Threshold; }

		virtual void setKnee(double value) { if (Knee != value) { Knee = value; setDirty(true); } }
		double getKnee() const { return Knee; }

		virtual void setRatio(double value) { if (Ratio != value) { Ratio = value; setDirty(true); } }
		double getRatio() const { return Ratio; }
		
		virtual void setMakeup(double value) { if (Makeup != value) { Makeup = value; setDirty(true); } }
		double getMakeup() const { return Makeup; }

		virtual void setMix(double value) { if (Mix != value) { Mix = value; setDirty(true); } }
		double getMix() const { return Mix; }

		// get/set Attributes
		virtual void setBackColor(CColor color) { if (BackColor != color) { BackColor = color; setDirty(true); } }
		CColor getBackColor() const { return BackColor; }

		virtual void setLineColor(CColor color) { if (LineColor != color) { LineColor = color; setDirty(true); } }
		CColor getLineColor() const { return LineColor; }


		// overrides
		void setDirty(bool state) override { CView::setDirty(state); };
		void draw(CDrawContext* _pContext) override {

			CDrawContext* pContext = _pContext;

			// draw border
			pContext->setLineWidth(1);
			pContext->setFillColor(BackColor);
			pContext->setFrameColor(VSTGUI::CColor(223, 233, 233, 255)); // black borders
			pContext->drawRect(getViewSize(), VSTGUI::kDrawFilledAndStroked);

			// draw db lines
			
			//double FREQ_LOG_MAX = log(MAX_FREQ / MIN_FREQ);
			//double DB_EQ_RANGE = 15.0;
			{
				VSTGUI::CRect r(getViewSize());
				auto width = r.getWidth();
				auto height = r.getHeight();

				pContext->setFrameColor(VSTGUI::CColor(255, 255, 255, 55));
				for (int x = MIN_dB; x < MAX_dB; x += 10) {
					VSTGUI::CCoord ver = width * -x * Inv_DB_R;
					const VSTGUI::CPoint _p1(r.left + ver, r.bottom);
					const VSTGUI::CPoint _p2(r.left + ver, r.top);
					pContext->drawLine(_p1, _p2);
				}

				for (int x = MIN_dB; x < MAX_dB; x += 10) {
					VSTGUI::CCoord ver = height * -x * Inv_DB_R;
					const VSTGUI::CPoint _p1(r.left, r.bottom - ver);
					const VSTGUI::CPoint _p2(r.right, r.bottom - ver);
					pContext->drawLine(_p1, _p2);
				}
			}


			// draw curve
			VSTGUI::CGraphicsPath* path = pContext->createGraphicsPath();
			if (path)
			{
				VSTGUI::CRect r(getViewSize());
				// VSTGUI::CCoord inset = 30;
				// r.inset(inset, 0);

				double thrh = getThrershold() * (0.0 - (-60.0)) + (-60.0);
				double knee = getKnee() * 20.0;
				double kneh = knee / 2.0;
				double inv_knee = 1.0 / knee;
				double ratio = getRatio() * (20.0 - (1.0)) + (1.0);
				double slope = 1.0 / ratio - 1.0;
				double makeup = getMakeup() * (12.0 - (-12.0)) + (-12.0);
				double width = r.getWidth();
				double inv_width = 1.0 / r.getWidth();
				double height = r.getHeight();

				path->beginSubpath(VSTGUI::CPoint(r.left - 1, r.bottom));
				for (int x = -1; x <= width + 1; x++) {
					double in = -DB_Range * (1.0 - x * inv_width); // -60 -> 0
					double overshoot = in - thrh;

					double y = 0.0;
					if (overshoot <= -kneh)
						y = 0.0;
					else if (overshoot > -kneh && overshoot <= kneh)
						y = 0.5 * slope * ((overshoot + kneh) * (overshoot + kneh) * inv_knee);
					else 
						y = slope * overshoot;

					double in_ = std::pow(10.0, (in) * (0.05)); // 0.001 -> 1
					double re_ = std::pow(10.0, (in + y + makeup) * (0.05));
					double m = ((1.0 - getMix()) * in_ + getMix() * re_); // 0.001 -> 1
					// double m = re_; // 0.001 -> 1

					double mm = (std::log10(m)) * (20.0); // -60 ~ 0
					double mmm = -mm * Inv_DB_R; // 1 ~ 0
					double scy = (1.0 - mmm) * height;
					path->addLine(VSTGUI::CPoint(r.left + x, r.bottom - scy));
				}
				path->addLine(VSTGUI::CPoint(r.right + 1, r.bottom + 1));
				path->addLine(VSTGUI::CPoint(r.left - 1, r.bottom + 1));
				path->closeSubpath();
				
				pContext->setFrameColor(LineColor);
				pContext->setDrawMode(VSTGUI::kAntiAliasing);
				pContext->setLineWidth(1.5);
				pContext->setLineStyle(VSTGUI::kLineSolid);
				pContext->drawGraphicsPath(path, VSTGUI::CDrawContext::kPathStroked);
				path->forget();
				
			}

			setDirty(false);
		};
		void setViewSize(const CRect& newSize, bool invalid = true) override
		{
			CControl::setViewSize(newSize, invalid);
		};
		bool sizeToFit() override {
			if (getDrawBackground())
			{
				CRect vs(getViewSize());
				vs.setWidth(getDrawBackground()->getWidth());
				vs.setHeight(getDrawBackground()->getHeight());
				setViewSize(vs);
				setMouseableArea(vs);
				return true;
			}
			return false;
		};
		/*
		void onIdle() override {
			if (getOldValue() != value)
				invalid();
		};
		*/

		CLASS_METHODS(CurveView, CControl)

	protected:
		~CurveView() noexcept override
		{
			//setOnBitmap(nullptr);
			//setOffBitmap(nullptr);
		};

		CColor		BackColor;
		CColor		LineColor;

		double Threshold = 0.0;
		double Knee = 0.0;
		double Ratio = 0.0;
		double Makeup = 0.0;
		double Mix = 1.0;

		const double MAX_dB = 0.0;
		const double MIN_dB = -60.0;
		const double DB_Range = MAX_dB - MIN_dB;
		const double Inv_DB_R = 1.0 / DB_Range;
	};


	//------------------------------------------------------------------------
	//  Metering reset container
	//------------------------------------------------------------------------
	class MeterViewContainer : public CViewContainer {
	public:
		MeterViewContainer(const CRect& size) : CViewContainer(size) {};
		void rr() {
			for (auto& child : getChildren())
			{
				// proc(child);
				// child.cast<CControl>()->setValue(child.cast<CControl>()->getDefaultValue());
				if (auto pp = child.cast<MyVuMeter>()) {
					//if (pp->getTag() == kInL)
						//pp->setValue(pp->getDefaultValue());
					pp->setValue(pp->getDefaultValue());
					pp->setDirty(true);
				}
				if (auto pp = child.cast<CParamDisplay>()) {
					pp->CParamDisplay::setValue(pp->getDefaultValue());
					pp->setDirty();
				}
			}
		}
		void onMouseDownEvent(MouseDownEvent& event) override {
			double test = 0.1;
			// forEachChild (Proc proc);
			for (auto& child : getChildren())
			{
				if (auto pp = child.cast<MeterViewContainer>()) {
					pp->rr();
				}
			}
			CViewContainer::onMouseDownEvent(event);
		};
		
	};
}

namespace yg331 {
	//------------------------------------------------------------------------
	// VuMeterController
	//------------------------------------------------------------------------
	template <typename ControllerType>
	class VuMeterController : public VSTGUI::IController, public VSTGUI::ViewListenerAdapter
	{
	public:
		VuMeterController(ControllerType* _mainController) : 
			mainController(_mainController), 
			inMeter(nullptr),
			outMeter(nullptr),
			grMeter(nullptr),
			vuMeterInL(nullptr),
			vuMeterInR(nullptr),
			vuMeterOutL(nullptr),
			vuMeterOutR(nullptr),
			vuMeterGR(nullptr)
		{
		}
		~VuMeterController() override
		{
			if (inMeter)     viewWillDelete(inMeter);
			if (outMeter)    viewWillDelete(outMeter);
			if (grMeter)     viewWillDelete(grMeter);
			if (vuMeterInL)  viewWillDelete(vuMeterInL);
			if (vuMeterInR)  viewWillDelete(vuMeterInR);
			if (vuMeterOutL) viewWillDelete(vuMeterOutL);
			if (vuMeterOutR) viewWillDelete(vuMeterOutR);
			if (vuMeterGR)   viewWillDelete(vuMeterGR);

			mainController->removeUIMessageController(this);
		}

		void setMeterValue(double val, int32_t tag)
		{
			if (tag == kIn) {
				if (!inMeter) return;
				inMeter->setValue(val);
				inMeter->setDirty(true);
			}
			if (tag == kOut) {
				if (!outMeter) return;
				outMeter->setValue(val);
				outMeter->setDirty(true);
			}
			if (tag == kGR) {
				if (!grMeter) return;
				grMeter->setValue(val);
				grMeter->setDirty(true);
			}
		}
		void setVuMeterValue(double inL, double inR, double outL, double outR, double GR)
		{
			if (vuMeterInL && vuMeterInR) {
				vuMeterInL->setValue(inL);
				vuMeterInR->setValue(inR);
				vuMeterInL->setDirty(true);
				vuMeterInR->setDirty(true);
			}
			else if (vuMeterOutL && vuMeterOutR) {
				vuMeterOutL->setValue(outL);
				vuMeterOutR->setValue(outR);
				vuMeterOutL->setDirty(true);
				vuMeterOutR->setDirty(true);
			}
			else if (vuMeterGR) {
				vuMeterGR->setValue(GR);
				vuMeterGR->setDirty(true);
			}
		}
		void setVuMeterValue(double val, int32_t tag)
		{
			if (tag == kInL) {
				if (!vuMeterInL) return;
				vuMeterInL->setValue(val);
				vuMeterInL->setDirty(true);
			}
			if (tag == kInR) {
				if (!vuMeterInR) return;
				vuMeterInR->setValue(val);
				vuMeterInR->setDirty(true);
			}
			if (tag == kOutL) {
				if (!vuMeterOutL) return;
				vuMeterOutL->setValue(val);
				vuMeterOutL->setDirty(true);
			}
			if (tag == kOutR) {
				if (!vuMeterOutR) return;
				vuMeterOutR->setValue(val);
				vuMeterOutR->setDirty(true);
			}
			if (tag == kGR) {
				if (!vuMeterGR) return;
				vuMeterGR->setValue(val);
				vuMeterGR->setDirty(true);
			}
		}

		

	private:
		using CControl = VSTGUI::CControl;
		using CView = VSTGUI::CView;
		using CParamDisplay = VSTGUI::CParamDisplay;
		using MyPD = VSTGUI::MyPD;
		using MyVuMeter = VSTGUI::MyVuMeter;
		using UTF8String = VSTGUI::UTF8String;
		using UIAttributes = VSTGUI::UIAttributes;
		using IUIDescription = VSTGUI::IUIDescription;

		//--- from IControlListener ----------------------
		void valueChanged (CControl* /*pControl*/) override {}
		void controlBeginEdit(CControl* /*pControl*/) override {}
		void controlEndEdit(CControl* pControl) override {}
		//--- is called when a view is created -----
		CView* verifyView(
			CView* view, 
			const UIAttributes& /*attributes*/,
			const IUIDescription* /*description*/) override
		{
			if (MyPD* control = dynamic_cast<MyPD*>(view); control)
			{
				if (control->getTag() == kIn) {
					inMeter = control;
					inMeter->registerViewListener(this);
					inMeter->CParamDisplay::setValue(inMeter->getDefaultValue());
				}
				if (control->getTag() == kOut) {
					outMeter = control;
					outMeter->registerViewListener(this);
					outMeter->CParamDisplay::setValue(outMeter->getDefaultValue());
				}
				if (control->getTag() == kGR) {
					grMeter = control;
					grMeter->registerViewListener(this);
					grMeter->CParamDisplay::setValue(grMeter->getDefaultValue());
				}
			}
			
			if (MyVuMeter* control = dynamic_cast<MyVuMeter*>(view); control) {
				if (control->getTag() == kInL) {
					vuMeterInL = control;
					vuMeterInL->registerViewListener(this);
					///vuMeterInL->setValue(0.0);
				}
				if (control->getTag() == kInR) {
					vuMeterInR = control;
					vuMeterInR->registerViewListener(this);
					//vuMeterInR->setValue(0.0);
				}
				if (control->getTag() == kOutL) {
					vuMeterOutL = control;
					vuMeterOutL->registerViewListener(this);
					//vuMeterOutL->setValue(0.0);
				}
				if (control->getTag() == kOutR) {
					vuMeterOutR = control;
					vuMeterOutR->registerViewListener(this);
					//vuMeterOutR->setValue(0.0);
				}
				if (control->getTag() == kGR) {
					vuMeterGR = control;
					vuMeterGR->registerViewListener(this);
					//vuMeterGR->setValue(0.0);
				}
			}
			
			return view;
		}
		//--- from IViewListenerAdapter ----------------------
		//--- is called when a view will be deleted: the editor is closed -----
		void viewWillDelete(CView* view) override
		{
			if (dynamic_cast<CParamDisplay*> (view) == inMeter && inMeter)
			{
				inMeter->unregisterViewListener(this);
				inMeter = nullptr;
			}
			if (dynamic_cast<CParamDisplay*> (view) == outMeter && outMeter)
			{
				outMeter->unregisterViewListener(this);
				outMeter = nullptr;
			}
			if (dynamic_cast<CParamDisplay*> (view) == grMeter && grMeter)
			{
				grMeter->unregisterViewListener(this);
				grMeter = nullptr;
			}
			
			if (dynamic_cast<MyVuMeter*>(view) == vuMeterInL && vuMeterInL) {
				vuMeterInL->unregisterViewListener(this);
				vuMeterInL = nullptr;
			}
			if (dynamic_cast<MyVuMeter*>(view) == vuMeterInR && vuMeterInR) {
				vuMeterInR->unregisterViewListener(this);
				vuMeterInR = nullptr;
			}
			if (dynamic_cast<MyVuMeter*>(view) == vuMeterOutL && vuMeterOutL) {
				vuMeterOutL->unregisterViewListener(this);
				vuMeterOutL = nullptr;
			}
			if (dynamic_cast<MyVuMeter*>(view) == vuMeterOutR && vuMeterOutR) {
				vuMeterOutR->unregisterViewListener(this);
				vuMeterOutR = nullptr;
			}
			if (dynamic_cast<MyVuMeter*>(view) == vuMeterGR && vuMeterGR) {
				vuMeterGR->unregisterViewListener(this);
				vuMeterGR = nullptr;
			}
		}

		ControllerType* mainController;
		CParamDisplay* inMeter;
		CParamDisplay* outMeter;
		CParamDisplay* grMeter;
		MyVuMeter* vuMeterInL;
		MyVuMeter* vuMeterInR;
		MyVuMeter* vuMeterOutL;
		MyVuMeter* vuMeterOutR;
		MyVuMeter* vuMeterGR;
	};



	//------------------------------------------------------------------------
	// CurveController
	//------------------------------------------------------------------------
	template <typename ControllerType>
	class CurveController 
		: public VSTGUI::ViewListenerAdapter
		, public Steinberg::FObject
		, public VSTGUI::DelegationController
	{
	public:
		CurveController(
			IController* baseController,
			ControllerType* _mainController, 
			Steinberg::Vst::Parameter* _paramThreshold, 
			Steinberg::Vst::Parameter* _paramKnee,
			Steinberg::Vst::Parameter* _paramRatio,
			Steinberg::Vst::Parameter* _paramMakeup,
			Steinberg::Vst::Parameter* _paramMix
		)
			: DelegationController(baseController)
			, mainController(_mainController)
			, paramThreshold(_paramThreshold)
			, paramKnee(_paramKnee)
			, paramRatio(_paramRatio)
			, paramMakeup(_paramMakeup)
			, paramMix(_paramMix)
			, curveView(nullptr) 
		{
			if (paramThreshold) paramThreshold->addDependent(this);
			if (paramKnee     ) paramKnee     ->addDependent(this);
			if (paramRatio    ) paramRatio    ->addDependent(this);
			if (paramMakeup   ) paramMakeup   ->addDependent(this);
			if (paramMix      ) paramMix      ->addDependent(this);
		}

		~CurveController() override
		{
			if (paramThreshold) paramThreshold->removeDependent(this);
			if (paramKnee     ) paramKnee     ->removeDependent(this);
			if (paramRatio    ) paramRatio    ->removeDependent(this);
			if (paramMakeup   ) paramMakeup   ->removeDependent(this);
			if (paramMix      ) paramMix      ->removeDependent(this);

			if (curveView) viewWillDelete(curveView);

			mainController->removeUICurveController(this);
		}

	private:
		using CControl = VSTGUI::CControl;
		using CView = VSTGUI::CView;
		using CurveView = VSTGUI::CurveView;
		using UTF8String = VSTGUI::UTF8String;
		using UIAttributes = VSTGUI::UIAttributes;
		using IUIDescription = VSTGUI::IUIDescription;

		void PLUGIN_API update(
			Steinberg::FUnknown* changedUnknown,
			Steinberg::int32 message)
		{
			if (curveView)
			{
				auto* p = Steinberg::FCast<Steinberg::Vst::Parameter>(changedUnknown);
				if (p 
					&& (p == paramThreshold) 
					|| (p == paramKnee) 
					|| (p == paramRatio) 
					|| (p == paramMakeup) 
					|| (p == paramMix))
				{
					if (message == kChanged)
					{
						if (p == paramThreshold && paramThreshold) curveView->setThreshold(p->getNormalized());
						if (p == paramKnee      && paramKnee     ) curveView->setKnee     (p->getNormalized());
						if (p == paramRatio     && paramRatio    ) curveView->setRatio    (p->getNormalized());
						if (p == paramMakeup    && paramMakeup   ) curveView->setMakeup   (p->getNormalized());
						if (p == paramMix       && paramMix      ) curveView->setMix      (p->getNormalized());
						curveView->invalid();
					}
					else if (message == kWillDestroy)
					{
						if (paramThreshold) paramThreshold->removeDependent(this);
						if (paramKnee     ) paramKnee     ->removeDependent(this);
						if (paramRatio    ) paramRatio    ->removeDependent(this);
						if (paramMakeup   ) paramMakeup   ->removeDependent(this);
						if (paramMix      ) paramMix      ->removeDependent(this);
						paramThreshold = nullptr;
						paramKnee      = nullptr;
						paramRatio     = nullptr;
						paramMakeup    = nullptr;
						paramMix       = nullptr;
					}
				}
			}
		}

		//--- from IControlListener ----------------------
		void valueChanged(CControl* pControl) override {
			/*
			if (pControl == curveView && curveView)
			{
				//float x, y;
				//CXYPad::calculateXY(pControl->getValue(), x, y);

				auto xId = paramThreshold->getInfo().id;
				if (mainController->setParamNormalized(xId, Threshold) == Steinberg::kResultTrue)
					mainController->performEdit(xId, mainController->getParamNormalized(xId));

				auto yId = paramKnee->getInfo().id;
				if (mainController->setParamNormalized(yId, y) == Steinberg::kResultTrue)
					mainController->performEdit(yId, mainController->getParamNormalized(yId));
			}
			*/
		}
		void controlBeginEdit(CControl* /*pControl*/) override {}
		void controlEndEdit(CControl* pControl) override {
			/*
			if (pControl == curveView && paramThreshold && paramKnee)
			{
				editController->endEdit(paramThreshold->getInfo().id);
				editController->endEdit(paramKnee->getInfo().id);
				editController->finishGroupEdit();
			}
			else
			{
				DelegationController::controlEndEdit(pControl);
			}
			*/
		}

		//--- is called when a view is created -----
		CView* verifyView(
			CView* view,
			const UIAttributes& /*attributes*/,
			const IUIDescription* /*description*/) override
		{
			if (CurveView* control = dynamic_cast<CurveView*>(view); control)
			{
				curveView = control;
				curveView->registerControlListener(this);
				//curveView->setValue(0.0);

				update(paramKnee,      kChanged);
				update(paramThreshold, kChanged);
				update(paramRatio,     kChanged);
				update(paramMakeup,    kChanged);
				update(paramMix,       kChanged);
			}
			return view;
		}
		//--- from IViewListenerAdapter ----------------------
		//--- is called when a view will be deleted: the editor is closed -----
		void viewWillDelete(CView* view) override
		{
			if (dynamic_cast<CurveView*> (view) == curveView && curveView)
			{
				curveView->unregisterControlListener(this);
				curveView = nullptr;
			}
		}

		ControllerType* mainController;
		Steinberg::Vst::Parameter* paramKnee;
		Steinberg::Vst::Parameter* paramThreshold;
		Steinberg::Vst::Parameter* paramRatio;
		Steinberg::Vst::Parameter* paramMakeup;
		Steinberg::Vst::Parameter* paramMix;
		CurveView* curveView;
	};

	
//------------------------------------------------------------------------
//  comp_Controller
//------------------------------------------------------------------------
class comp_Controller 
	: public Steinberg::Vst::EditControllerEx1
	, public VSTGUI::VST3EditorDelegate
	, public Steinberg::Vst::IDataExchangeReceiver
{
public:
	using UIMessageController = VuMeterController<comp_Controller>;
	using UICurveController = CurveController<comp_Controller>;

	comp_Controller () = default;
	~comp_Controller () SMTG_OVERRIDE = default;

    // Create function
	static Steinberg::FUnknown* createInstance (void* /*context*/)
	{
		return (Steinberg::Vst::IEditController*)new comp_Controller;
	}

	// IPluginBase
	Steinberg::tresult PLUGIN_API initialize (Steinberg::FUnknown* context) SMTG_OVERRIDE;
	Steinberg::tresult PLUGIN_API terminate () SMTG_OVERRIDE;

	// EditController
	Steinberg::tresult PLUGIN_API setComponentState (Steinberg::IBStream* state) SMTG_OVERRIDE;
	Steinberg::IPlugView* PLUGIN_API createView (Steinberg::FIDString name) SMTG_OVERRIDE;
	Steinberg::tresult PLUGIN_API setState (Steinberg::IBStream* state) SMTG_OVERRIDE;
	Steinberg::tresult PLUGIN_API getState (Steinberg::IBStream* state) SMTG_OVERRIDE;
	Steinberg::tresult PLUGIN_API setParamNormalized (Steinberg::Vst::ParamID tag,
                                                      Steinberg::Vst::ParamValue value) SMTG_OVERRIDE;
	Steinberg::tresult PLUGIN_API getParamStringByValue (Steinberg::Vst::ParamID tag,
                                                         Steinberg::Vst::ParamValue valueNormalized,
                                                         Steinberg::Vst::String128 string) SMTG_OVERRIDE;
	Steinberg::tresult PLUGIN_API getParamValueByString (Steinberg::Vst::ParamID tag,
                                                         Steinberg::Vst::TChar* string,
                                                         Steinberg::Vst::ParamValue& valueNormalized) SMTG_OVERRIDE;
	/*
	VSTGUI::CView* createCustomView(VSTGUI::UTF8StringPtr name, const VSTGUI::UIAttributes& attributes,
		const VSTGUI::IUIDescription* description, VSTGUI::VST3Editor* editor) override
	{
		if (VSTGUI::UTF8StringView(name) == "EQCurveView")
		{
			VSTGUI::CRect size(VSTGUI::CPoint(30, 10), VSTGUI::CPoint(620, 150));
			return VSTGUI::EQCurveView(size, this);
		}
		return nullptr;
	}
	*/


	VSTGUI::CView* createCustomView(VSTGUI::UTF8StringPtr name, const VSTGUI::UIAttributes& attributes,
		const VSTGUI::IUIDescription* description, VSTGUI::VST3Editor* editor) override
	{
		return nullptr;
	};

	VSTGUI::CView* verifyView(
		VSTGUI::CView* view,
		const VSTGUI::UIAttributes& attributes,
		const VSTGUI::IUIDescription* description)
	{
		return view;
	};
	
	VSTGUI::IController* createSubController(
		VSTGUI::UTF8StringPtr name, 
		const VSTGUI::IUIDescription* description,
		VSTGUI::VST3Editor* editor) override
	{
		if (VSTGUI::UTF8StringView(name) == "MessageController")
		{
			auto* controller = new UIMessageController (this);
			addUIMessageController(controller);
			return controller;
		}
		if (VSTGUI::UTF8StringView(name) == "CurveController")
		{
			Steinberg::Vst::Parameter* ThresholdParam = getParameterObject(kParamThreshold);
			Steinberg::Vst::Parameter* KneeParam = getParameterObject(kParamKnee);
			Steinberg::Vst::Parameter* RatioParam = getParameterObject(kParamRatio);
			Steinberg::Vst::Parameter* MakeupParam = getParameterObject(kParamMakeup);
			Steinberg::Vst::Parameter* MixParam = getParameterObject(kParamMix);
			auto* controller = new UICurveController(editor, this, ThresholdParam, KneeParam, RatioParam, MakeupParam, MixParam);
			addUICurveController(controller);
			return controller;
		}
		return nullptr;
	};

	// EditController
	tresult PLUGIN_API notify(Vst::IMessage* message) override;

	void PLUGIN_API update(Steinberg::FUnknown* changedUnknown, Steinberg::int32 message) SMTG_OVERRIDE;
	//void editorAttached(Steinberg::Vst::EditorView* editor) SMTG_OVERRIDE;
	//void editorRemoved(Steinberg::Vst::EditorView* editor) SMTG_OVERRIDE;


	// IDataExchangeReceiver
	void PLUGIN_API queueOpened (Vst::DataExchangeUserContextID userContextID, uint32 blockSize,
	                             TBool& dispatchOnBackgroundThread) override;
	void PLUGIN_API queueClosed (Vst::DataExchangeUserContextID userContextID) override;
	void PLUGIN_API onDataExchangeBlocksReceived (Vst::DataExchangeUserContextID userContextID,
	                                              uint32 numBlocks, Vst::DataExchangeBlock* blocks,
	                                              TBool onBackgroundThread) override;

	//---Internal functions-------
	void addUICurveController(UICurveController* controller) {
		curveControllers.push_back(controller);
	};
	void removeUICurveController(UICurveController* controller) {
		auto it = std::find(curveControllers.begin(), curveControllers.end(), controller);
		if (it != curveControllers.end())
			curveControllers.erase(it);
	};

	void addUIMessageController(UIMessageController* controller) {
		vuMeterControllers.push_back(controller);
	};
	void removeUIMessageController(UIMessageController* controller) {
		auto it = std::find(vuMeterControllers.begin(), vuMeterControllers.end(), controller);
		if (it != vuMeterControllers.end())
			vuMeterControllers.erase(it);
	};

	void setDefaultMessageText(Vst::String128 text) {
		Steinberg::String tmp(text);
		tmp.copyTo16(defaultMessageText, 0, 127);
	};
	Vst::TChar* getDefaultMessageText() {
		return defaultMessageText;
	};

	double* getL() { return &InMeter[0]; };

 	//---Interface---------
	DEFINE_INTERFACES
		// Here you can add more supported VST3 interfaces
		DEF_INTERFACE(IDataExchangeReceiver)
	END_DEFINE_INTERFACES (EditController)
    DELEGATE_REFCOUNT (EditController)
//------------------------------------------------------------------------
protected:
	typedef std::vector<Steinberg::Vst::EditorView*> EditorVector;
	EditorVector editors;
	VSTGUI::VST3Editor* eded = nullptr;

	struct ZoomFactor {
		const Steinberg::tchar* title;
		double factor;

		ZoomFactor(const Steinberg::tchar* title, double factor) : title(title), factor(factor) {}
	};
	typedef std::vector<ZoomFactor> ZoomFactorVector;
	ZoomFactorVector zoomFactors;
	

	using UIMessageControllerList = std::vector<UIMessageController*>;
	UIMessageControllerList vuMeterControllers;
	using UICurveControllerList = std::vector<UICurveController*>;
	UICurveControllerList curveControllers;

	Vst::String128 defaultMessageText;

	double InMeter[2] = { 0.0, };
	double OutMeter[2] = { 0.0, };
	double GRMeter = { 0.0 };
	
	Vst::DataExchangeReceiverHandler dataExchange{ this };
};

//------------------------------------------------------------------------
} // namespace yg331
