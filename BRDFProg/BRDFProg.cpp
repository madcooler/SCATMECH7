//************************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: BRDFProg.cpp
//**
//** Thomas A. Germer
//** Optical Technology Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//** Version: 7.00 (January 2015)
//**
//************************************************************************************
#include "brdf.h"
#include "facet.h"
#include "twoface.h"
#include "allrough.h"
#include <iostream>
#include <fstream>

using namespace std;        // Use unqualified names for Standard C++ library
using namespace SCATMECH;   // Use unqualified names for SCATMECH library

int main(int argv,char** argc)
{
    try {

		// Query user for scattering angles and ranges...
		/* double thetai = AskUser("Incident Angle ",45.)*deg;

		double thetas_1 = AskUser("Initial Scattering Angle ",45.)*deg;
		double thetas_2 = AskUser("Final Scattering Angle ",45.)*deg;
		double thetas_3 = AskUser("Step Scattering Angle ",1.)*deg;

		double phis_1 = AskUser("Initial Azimuthal Angle ",0.)*deg;
		double phis_2 = AskUser("Final Azimuthal Angle ",180.)*deg;
		double phis_3 = AskUser("Step Azimuthal Angle ",2.)*deg;*/
		//const double 

		char filename [] = "Aluminium_RMS=0.1_635nm.txt";
		char filenameDielectric [] = "BRDF_Data.txt";
		ofstream ofs ( filenameDielectric );

		int thetai = 42;

		int thetas_1 = 32;
		int thetas_2 = 57;
		int thetas_3 = 5;

		int phis_1 = 0;
		int phis_2 = 0;
		int phis_3 = 1;

		// Get an instance of BRDF_Model...
		//BRDF_Model_Ptr model = Get_Model_Ptr();
		//model->AskUser();

		/*Roughness_BRDF_Model_Ptr model="Microroughness_BRDF_Model";

		model->set_lambda(0.635);
		optical_constant oc(1.37,7.62);
		model->set_substrate(oc);

		PSD_Function_Ptr psd="Gaussian_PSD_Function";
		psd->set_parameter("sigma",10);
		psd->set_parameter("length",0.5);
		model->set_psd(psd);*/

		//optical_constant layer1(1.3,0);
		//dielectric_function df2(layer2);
		//dielectric_function df1(layer1);
		//model->set_parameter("film",df1);



		//Roughness_BRDF_Model_Ptr model="Uncorrelated_Roughness_Stack_BRDF_Model"; 
		Facet_BRDF_Model_Ptr model = "Shadowed_Facet_BRDF_Model";
		//Facet_BRDF_Model_Ptr model="Subsurface_Facet_BRDF_Model";

		model->set_lambda ( 0.635 );
		optical_constant oc ( 0.23 , 3.46 );
		model->set_substrate ( oc );

		/*PSD_Function_Ptr psd="Gaussian_PSD_Function";
		psd->set_parameter("sigma",10);
		psd->set_parameter("length",0.5);
		model->set_psd(psd);*/



		Slope_Distribution_Function_Ptr sdf = "Gaussian_Slope_Distribution_Function";
		sdf->set_parameter ( "s" , 0.3 );
		model->set_sdf ( sdf );

		/*optical_constant top(1.0,0);
		model->set_parameter("overcoat",top);
		*/
		//optical_constant layer2(1.9,0);
		//optical_constant layer1(1.3,0);
		//dielectric_function df2(layer2);
		//dielectric_function df1(layer1);
		//dielectric_stack ds;
		//ds.grow(df2,0.1);
		//ds.grow(df1,0.1);
		//model->set_parameter("stack",ds);


		/*optical_constant layer2(6.5,0);
		optical_constant layer1(1.3,0);
		dielectric_function df2(layer2);
		dielectric_function df1(layer1);
		dielectric_stack ds;
		ds.grow(df2,1.5);
		ds.grow(df1,0.1);
		model->set_parameter("stack",ds);*/


		/*Facet_BRDF_Model_Ptr model="Shadowed_Facet_BRDF_Model";
		//Facet_BRDF_Model_Ptr model="Subsurface_Facet_BRDF_Model";
		Slope_Distribution_Function_Ptr sdf="Gaussian_Slope_Distribution_Function";
		//Inheritance sdfinheritance=sdf->get_inheritance();
		sdf->set_parameter("s",0.2);
		//sdf->f(0.2);

		// Query user for model parameters...
		//model->AskUser();
		model->set_sdf(sdf);
		//model->lambda=0.4;
		//model->set_parameter()
		optical_constant oc(1.37,7.62);
		//optical_constant oc(1.5,0);
		model->set_lambda(0.635);
		model->set_substrate(oc);
		*/
		// Loop through scattering geometries...
		////for (double thetai=0;thetai<=90*deg;thetai+=thetas_3) {
		for (int thetas = thetas_1; thetas <= thetas_2; thetas += thetas_3) {
			for (int phis = phis_1; phis <= phis_2; phis += phis_3) {

				// Get the Mueller matrix for scattering...
				MuellerMatrix m = model->Mueller ( thetai*deg , thetas*deg , phis*deg , 0 );
				//float t=thetas*pi/180;
				//float costhetas=cos(t);
				// Calculate the Stokes vector for p-polarized incident light...
				//MyStokesVector s=m*MyStokesVectorUnitUnpol();
				StokesVector incidentray = StokesVector ( 1 , 0 , 1 , 0 );
				StokesVector s = m * incidentray;
				//MyStokesVector s=m*MyStokesVectorUnitLCP();

				// Print out various light scattering parameters...
				/*std::cout << thetas/deg << tab   // Scattering angle (theta)
				<< phis/deg << tab     // Scattering angle (phi)
				//<< s.eta()/deg << tab  // Principal angle of polarization
				//<< s.DOLP() << tab     // Degree of linear polarization
				//<< s.DOCP() << tab     // Degree of circular polarization
				<< s.DOP() << tab      // Degree of polarization
				<< s.I() <<tab
				<<s.Q() << tab
				<<s.U() << tab
				<<s.V() << tab
				<<s.intensity()<<tab
				<< std::endl;      // The intensity (BRDF)*/



				ofs << thetas << " "   // Scattering angle (theta)
					<< phis << " "      // Scattering angle (phi)
										//<< s.eta()/deg << tab  // Principal angle of polarization
										//<< s.DOLP() << tab     // Degree of linear polarization
										//<< s.DOCP() << tab     // Degree of circular polarization

					<< s.I () << " "
					<< s.Q () << " "
					<< s.U () << " "
					<< s.V () << " " //<< s.DOP() <<" "      // Degree of polarization
									 //<<s.intensity()<<" " 
					<< s.DOP () << " "
					<< s.DOLP () << " "      // Degree of linear polarization
					<< s.DOCP () << " "      // Degree of circular polarization
					<< s.eta () / deg << " "   // Principal angle of polarization
					<< std::endl;      // The intensity (BRDF)

									   //ofs //<< thetas/deg <<" "   // Scattering angle (theta)
									   //<< phis /deg<<" "      // Scattering angle (phi)
									   //<< s.eta()/deg << tab  // Principal angle of polarization
									   //<< s.DOLP() << tab     // Degree of linear polarization
									   //<< s.DOCP() << tab     // Degree of circular polarization

									   //<< s.I() <<" " 
									   //<<s.Q() <<" "      
									   //<<s.U() <<" "  
									   //<<s.V() <<" " //<< s.DOP() <<" "      // Degree of polarization
									   //<<s.intensity()<<" " 
									   //<< s.DOP() <<" "
									   //<< s.DOLP() <<" "      // Degree of linear polarization
									   //<< s.DOCP() <<" "      // Degree of circular polarization
									   //<< s.eta()/deg <<" "   // Principal angle of polarization
									   //<< std::endl;      // The intensity (BRDF)

			}
		}

		return 0;

    } catch (exception& e) {

        cerr << e.what() << endl;

    }

    return 0;
}
