// DGP 2019 Project
// ShapeUp and Projective Dynamics
// Author: Christopher Brandt

// Modified widgets used in the ProjDynAPI

#pragma once

#include "projdyn_types.h"
#include "projdyn.h"
#include "projdyn_tetgen.h"
#include <nanogui/slider.h>
#include "viewer.h"

// A specialized constraint-weight slider widget, inheriting from the NanoGUI slider,
// which highlights the associated constraint groups on a mouse-over-event, and
// changes the weight of the constraint (with some exponential scaling).
// Comes with a text-box that shows the value.
class ConstraintSlider : public Slider {
public:
    ConstraintSlider(Widget* parent, Viewer* viewer, ProjDyn::Index numVerts, ProjDyn::ConstraintGroupPtr constraint)
        : Slider(parent) {
        m_constraint = constraint;
        m_viewer = viewer;
        m_numVerts = numVerts;

        if (constraint->weight > 0) {
            setValue(std::log(constraint->weight) / std::log(10000) + 1);
        }
        else {
            setValue(0);
        }
        setRange(std::pair<float, float>(0, 2));
        setFixedWidth(80);

        // Add a textbox and set defaults
        m_textBox = new TextBox(parent);
        m_textBox->setFixedSize(Vector2i(80, 25));
        m_textBox->setValue(ProjDyn::floatToString(value()));

        setCallback([this](float v) {
            float vv = std::pow(10000, v - 1);
            m_textBox->setValue(ProjDyn::floatToString(v));
            m_constraint->weight = vv;
        });

        m_textBox->setEditable(true);
        m_textBox->setCallback([this](const std::string& val) -> bool {
            float v = std::stof(val);
            float vv = std::pow(10000, v - 1);
            m_constraint->weight = vv;
            this->setValue(v);
            this->finalCallback()(v);
            return true;
        });

    }

    virtual bool mouseEnterEvent(const Vector2i& p, bool enter) override {
        if (enter) {
            m_entered = true;
            m_storedStatus = m_viewer->getVertexStatus();
            Eigen::RowVectorXi tempStatus;
            tempStatus.setZero(m_numVerts);
            for (auto c : m_constraint->constraints) {
                for (ProjDyn::Index ind : c->getIndices()) {
                    if (ind < m_numVerts) tempStatus(0, ind) = 1;
                }
            }
            m_viewer->updateVertexStatus(tempStatus);
            return true;
        }
        else if (m_entered) {
            m_entered = false;
            m_viewer->updateVertexStatus(m_storedStatus);
        }
        return false;
    }

    // This needs to be overridden since the event of leaving the area of the slider with the mouse
    // while clicking does not trigger the mouseEnterEvent, i.e. the vertex status is never
    // restored
    virtual bool mouseButtonEvent(const Vector2i& p, int button, bool down, int modifiers) override
    {
        Slider::mouseButtonEvent(p, button, down, modifiers);
        // When releasing the button and we "entered" this constraint before, restore the original
        // vertex status.
        if (!down && m_entered) {
            m_entered = false;
            m_viewer->updateVertexStatus(m_storedStatus);
        }
        return mEnabled;
    }

private:
    bool m_entered = false;
    ProjDyn::ConstraintGroupPtr m_constraint;
    Viewer* m_viewer;
    ProjDyn::Index m_numVerts;
    Eigen::RowVectorXi m_storedStatus;
    TextBox* m_textBox;
};