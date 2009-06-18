//
//  HelloAccelerometerController.h
//  HelloAccelerometer
//
//  Created by turner on 5/6/09.
//  Copyright Douglass Turner Consulting 2009. All rights reserved.
//


#import <UIKit/UIKit.h>
#import "VectorMatrix.h"

@interface HelloAccelerometerController : NSObject <UIAccelerometerDelegate> {
	
	// Raw acceleration straight from the device
    IBOutlet UILabel		*acc_raw_x_label;
    IBOutlet UILabel		*acc_raw_y_label;
    IBOutlet UILabel		*acc_raw_z_label;
	M3DVector3f				acc_raw;
	
    IBOutlet UILabel		*acc_raw_filtered_x_label;
    IBOutlet UILabel		*acc_raw_filtered_y_label;
    IBOutlet UILabel		*acc_raw_filtered_z_label;
	M3DVector3f				acc_raw_filtered;
	M3DVector3f				acc_raw_spherical_filtered;
	
	// Physical "G" value
    IBOutlet UILabel		*g_x_label;
    IBOutlet UILabel		*g_y_label;
    IBOutlet UILabel		*g_z_label;
	M3DVector3f				g;
	M3DVector3f				gPast;
	
	// Robert Paul notation for matrices
	//	nx ox ax px
	//	ny oy ay py
	//	nz oz az pz
    IBOutlet UILabel		*nx_label;
    IBOutlet UILabel		*ny_label;
    IBOutlet UILabel		*nz_label;
	M3DVector3f				x;
	
    IBOutlet UILabel		*ox_label;
    IBOutlet UILabel		*oy_label;
    IBOutlet UILabel		*oz_label;
	M3DVector3f				y;
	
    IBOutlet UILabel		*ax_label;
    IBOutlet UILabel		*ay_label;
    IBOutlet UILabel		*az_label;
	M3DVector3f				z;
	
	M3DMatrix44f			_cameraTransform;
	M3DMatrix44f			_openGLCameraInverseTransform;
	
	// Spherical
    IBOutlet UILabel		*rho_label;
    IBOutlet UILabel		*theta_label;
    IBOutlet UILabel		*phi_label;
	M3DVector3f				gSpherical;
	M3DVector3f				gSphericalPast;
	
	// Spherical deltas
    IBOutlet UILabel		*delta_rho_label;
    IBOutlet UILabel		*delta_theta_label;
    IBOutlet UILabel		*delta_phi_label;
	M3DVector3f				deltaSpherical;
	
	// Instantaneous acceleration
    IBOutlet UILabel		*acc_instantaneous_x_label;
    IBOutlet UILabel		*acc_instantaneous_y_label;
    IBOutlet UILabel		*acc_instantaneous_z_label;
	M3DVector3f				acc_instantaneous;
	
	NSMutableArray* _x_samples;
	NSMutableArray* _y_samples;
	NSMutableArray* _z_samples;
	
	NSMutableArray* _rho_samples;
	NSMutableArray* _theta_samples;
	NSMutableArray* _phi_samples;
}

- (void)initializeAccelerometer;

@end
