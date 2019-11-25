// DGP 2019 Project
// ShapeUp and Projective Dynamics
// Author: Christopher Brandt

// Program entry point for the Projective Dynamics simulator. 
// Simply sets up the simulation using the ProjDynAPI and the viewer.

#include "viewer.h"
#include "projdyn_api.h"
#include <thread>
#include <chrono>

ProjDynAPI* pd_api_ptr = nullptr;

// Some call back functions for handling mesh grab events and loading new meshes
// from the viewer.
void projdyn_grab(const std::vector<ProjDyn::Index>& grabbedVerts, const std::vector<Vector3f>& grabPos) {
    if (pd_api_ptr) {
        pd_api_ptr->grab(grabbedVerts, grabPos);
    }
}
void projdyn_release_grab() {
    if (pd_api_ptr) {
        pd_api_ptr->releaseGrab();
    }
}
bool projdyn_preload(Viewer* viewer) {
    if (pd_api_ptr) {
        pd_api_ptr->stop();
        return true;
    }
    else {
        return false;
    }
}
bool projdyn_setmesh(Viewer* viewer) {
    if (pd_api_ptr) {
        return pd_api_ptr->setMesh(false);
    }
    else {
        return false;
    }
}

// Entry point for the program
int main(int argc, char** argv) {
	try {
		// Initialize the GUI engine
		nanogui::init();

		// Create a new mesh viewer app, which will add a screen to nanogui
		// The callback function is triggered when loading a new mesh and (re-)initializes the
		// projective dynamics simulator with the new vertices, faces and tetrahedrons
		nanogui::ref<Viewer> app = new Viewer("Projective Dynamics", nullptr, projdyn_setmesh, projdyn_preload);
		app->setGrabCallbacks(projdyn_grab, projdyn_release_grab);

        // Initialize the projective dynamics GUI
        ProjDynAPI pd_api(app);
        pd_api_ptr = &pd_api;
        pd_api.setGravity(9.81);
		pd_api.initSimulationGUI(true);

		app->drawAll();
		app->setVisible(true);

		// Start the main loop: keeps calling drawContents() of the viewer, until the window is closed
		nanogui::mainloop(-1);

		// Clean-up
		pd_api.stop();
		nanogui::shutdown();

		// Error handling:
	}
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
