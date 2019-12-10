#ifndef PROJDYN_TEMPERATURE_H
#define PROJDYN_TEMPERATURE_H

// DGP 2019 Project
// ShapeUp and Projective Dynamics
// Authors: Shad Durussel, Amandine Evard, Michele Vidulis

// This file contains the Temperature class, which implements the
// temperature models used to simulate temperature changes on the mesh.

#include "projdyn.h"
#include "projdyn_types.h"
#include "projdyn_common.h"

namespace ProjDyn {
	
	//enum TemperatureModel { none, uniform, linear, diffusion };

	//-------------------------------------------------------------
	// Base abstract class
	//-------------------------------------------------------------
	class Temperature {
	public:
		Temperature(Simulator* simulator, TemperatureModel model) :
			m_temperature_model(model),
			m_temperatures(0),
			m_simulator(simulator)
		{}

    /* Method called at every time step by step().
    Updates temperature depending on the temperature_model chosen */
    virtual void updateTemperature() = 0;

    /* Returns the coefficient controlling the temperature behavior
    in the model chosen. Its name differs in every derived class */
    virtual const Scalar getTemperatureCoefficient() = 0;

    const TemperatureModel getTemperatureModel() const {
		return m_temperature_model;
	}

    const Vector& getTemperatures() const {
		return m_temperatures;
	}

	protected:
		TemperatureModel m_temperature_model;
		Vector m_temperatures;
		Simulator* m_simulator;
	};

	//-------------------------------------------------------------
	// Dervied classes implementing different temperature models
	//-------------------------------------------------------------

	/* Temperature scales uniformly on the whole mesh. 
	The value can be controlled from the GUI */
	class TemperatureUniform : public Temperature {
	public:
		TemperatureUniform(Simulator* simulator):
		  Temperature(simulator, TemperatureModel::uniform) {}

		const Scalar getTemperatureCoefficient() override{
		  return m_uniform_temperature;
		}

		void updateTemperature() override {
		  Vector newTemp;
				newTemp.resize(m_simulator->getNumVerts());
				newTemp.setOnes();
				newTemp *= m_uniform_temperature;
				m_temperatures = newTemp;
		}

	protected:
		Scalar m_uniform_temperature;
	};

	/* Temperature scales linearly with height.
	The floor temperature can be controlled from the GUI,
	the reference value far from the floor is 0 (TODO: generalize)*/
	class TemperatureLinear : public Temperature {
	public:
		TemperatureLinear(Simulator* simulator):
		  Temperature(simulator, TemperatureModel::linear) {}

		const Scalar getTemperatureCoefficient() override{
		  return m_linear_bottom_temperature;
		}

		void updateTemperature() override {
			// TODO: improve implementation (add also m_linear_top_temperature control and get rid of meaningless coef)
    		Scalar coef = 20;
    		Scalar slope = (m_linear_bottom_temperature - m_linear_top_temperature) / m_simulator->getFloorHeight() * coef;
    		Scalar intercept = m_linear_bottom_temperature - slope * m_simulator->getFloorHeight();
    		Vector interceptTemp;
    		interceptTemp.resize(m_simulator->getNumVerts());
    		interceptTemp.setOnes();
    		interceptTemp *= intercept;
    		m_temperatures = interceptTemp +  slope * m_simulator->getPositions().col(1);

    		for (int i = 0; i < m_simulator->getNumVerts(); i++) {
    			if (m_temperatures[i] < 0) {
    				m_temperatures[i] = 0;
    			}
    			else if (m_temperatures[i] > m_linear_bottom_temperature) {
    				m_temperatures[i] = m_linear_bottom_temperature;
    			}
    		}
		}

	protected:
		Scalar m_linear_top_temperature = 0; // with linear increase the reference value is always 0 (TODO: generalize)
		Scalar m_linear_bottom_temperature;
	};


	class TemperatureDiffusion : public Temperature {
	public:
		TemperatureDiffusion(Simulator* simulator):
			Temperature(simulator, TemperatureModel::diffusion) {}

		const Scalar getTemperatureCoefficient() override{
			return m_diffusion_coefficient;
		}

		void updateTemperature() override {
		Vector groundHeight;
			groundHeight.resize(m_simulator->getNumVerts());
			groundHeight.setOnes();
			groundHeight *= m_simulator->getFloorHeight();

			Positions positions = m_simulator->getPositions();

			for (int i = 0; i < m_simulator->getNumVerts(); i++) {
				if (positions.col(1)[i] - groundHeight[i] < 0.00001) {
					m_temperatures[i] = m_floor_temperature;
				}
			}

			Vector new_temp;
			new_temp.resize(m_simulator->getNumVerts());
			new_temp.setZero();
			Scalar time_step = m_simulator->getTimeStep();
			for (int i = 0; i < m_simulator->getNumVerts(); i++) {
				Scalar ti = m_temperatures[i];
				std::vector<int> neig = m_simulator->getNeighbors().at(i);
				int nb_neighbors = neig.size();
				//std::cout << "nb_neighbors i = " << nb_neighbors <<"\n";
				Scalar new_temp_i = 0;
				for(int j: neig){
					Scalar tj = m_temperatures[j];
					new_temp_i += 1.0/nb_neighbors * (tj-ti)*m_diffusion_coefficient*time_step/((positions.row(i)- positions.row(j)).norm());
				}
				new_temp[i] = ti+new_temp_i;
				//std::cout << "new temp i = " << new_temp_i <<"\n";
				if (new_temp[i] > m_floor_temperature){new_temp[i] = m_floor_temperature;}
				if(new_temp[i] <0 ){new_temp[i] = 0;}
			}
			m_temperatures = new_temp;
		}

	protected:
		Scalar m_diffusion_coefficient;
		Scalar m_floor_temperature = 100; // TODO: generalize
	};
}

#endif // PROJDYN_TEMPERATURE_H