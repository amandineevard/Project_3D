// DGP 2019 Project
// ShapeUp and Projective Dynamics
// Author: Christopher Brandt

// Program entry point for ShapeUp-style shape deformation. 
// Simply sets up the simulator (local-global algorithm) using the 
// ShapeUpAPI and the viewer.

#pragma once

#include "viewer.h"
#include "projdyn_api.h"
#include "shapeup_api.h"
#include <thread>
#include <chrono>

ShapeUpAPI* shapeup_ptr = nullptr;

// The following are some global functions that handle the call-backs
// from the viewer for changes of the mesh and mouse events.
// They are re-routed to the current ShapeUpAPI object.
bool shapeup_set_mesh(Viewer* viewer) {
    if (shapeup_ptr) return shapeup_ptr->setMesh();
    return false;
}
void shapeup_grab(const std::vector<ProjDyn::Index>& grabbedVerts, const std::vector<Vector3f>& grabPos) {
    if (shapeup_ptr) shapeup_ptr->mouseGrabEvent(grabbedVerts, grabPos);
    return;
}
void shapeup_release_grab() {
    if (shapeup_ptr) shapeup_ptr->mouseGrabRelease();
}

// Entry point for the program
int main(int argc, char** argv) {
    try {
        // Initialize the GUI engine
		nanogui::init();

		// Create a new mesh viewer app, which will add a screen to nanogui
		// The callback function is triggered when loading a new mesh and initializes
        // the constraints used by ShapeUp
		nanogui::ref<Viewer> app = new Viewer("ShapeUp", nullptr, shapeup_set_mesh);

        // Initialize the ShapeUp GUI
        ShapeUpAPI shapeup_api(app);
        shapeup_ptr = &shapeup_api;

        // The following callback functions handle the response when alt-left-clicking and dragging
        // vertices in the viewer.
        app->setGrabCallbacks(shapeup_grab, shapeup_release_grab);

        // Add GUI controls for ShapeUp
        shapeup_api.initShapeUpGUI();

		app->drawAll();
		app->setVisible(true);

		// Start the main loop: keeps calling drawContents() of the viewer, until the window is closed
		nanogui::mainloop(-1);

		// Clean-up
		nanogui::shutdown();

	}
    // Error handling:
    catch (const std::runtime_error& e) {
		std::string error_msg = std::string("Caught a fatal error: ") + std::string(e.what());
#if defined(_WIN32)
		MessageBoxA(nullptr, error_msg.c_str(), NULL, MB_ICONERROR | MB_OK);
#else
		std::cerr << error_msg << std::endl;
#endif
		return -1;
	}

	return 0;
}