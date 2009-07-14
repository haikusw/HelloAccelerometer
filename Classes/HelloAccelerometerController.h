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

	NSMutableArray			*_x_samples;
	NSMutableArray			*_y_samples;
	NSMutableArray			*_z_samples;

	M3DVector3f				acc_instantaneous;

	M3DVector3f				g;
	M3DVector3f				gSpherical;
	
	M3DVector3f				gPast;

	M3DVector3f				exe;
	M3DVector3f				wye;
	M3DVector3f				zee;
	
	M3DMatrix44f			cameraTransform;
	M3DMatrix44f			modelingTransform;
	M3DMatrix44f			openGLModelViewTransform;

	
	// Robert Paul notation for matrices
	//	nx ox ax px
	//	ny oy ay py
	//	nz oz az pz
	
	// The camera transform.
    IBOutlet UILabel		*nx_label;
    IBOutlet UILabel		*ny_label;
    IBOutlet UILabel		*nz_label;
	
    IBOutlet UILabel		*ox_label;
    IBOutlet UILabel		*oy_label;
    IBOutlet UILabel		*oz_label;
	
    IBOutlet UILabel		*ax_label;
    IBOutlet UILabel		*ay_label;
    IBOutlet UILabel		*az_label;
	
    IBOutlet UILabel		*px_label;
    IBOutlet UILabel		*py_label;
    IBOutlet UILabel		*pz_label;

	
	// The modeling transform
    IBOutlet UILabel		*nx_model_label;
    IBOutlet UILabel		*ny_model_label;
    IBOutlet UILabel		*nz_model_label;
	
    IBOutlet UILabel		*ox_model_label;
    IBOutlet UILabel		*oy_model_label;
    IBOutlet UILabel		*oz_model_label;
	
    IBOutlet UILabel		*ax_model_label;
    IBOutlet UILabel		*ay_model_label;
    IBOutlet UILabel		*az_model_label;
	
    IBOutlet UILabel		*px_model_label;
    IBOutlet UILabel		*py_model_label;
    IBOutlet UILabel		*pz_model_label;
	
	
	
	
	
	
	
	
	
	// "G" Vector
    IBOutlet UILabel		*g_x_label;
    IBOutlet UILabel		*g_y_label;
    IBOutlet UILabel		*g_z_label;
	
	// "G" Vector in spherical coordinates
    IBOutlet UILabel		*rho_label;
    IBOutlet UILabel		*theta_label;
    IBOutlet UILabel		*phi_label;
	
	// Instantaneous acceleration
    IBOutlet UILabel		*acc_instantaneous_x_label;
    IBOutlet UILabel		*acc_instantaneous_y_label;
    IBOutlet UILabel		*acc_instantaneous_z_label;
	
}

+(NSString*)signPrint:(float)f;

+ (NSString*)signedAxisName:(M3DVector3f)vector axis:(NSUInteger)axis;

+(float)sign:(float)f;

+ (void)dominantAxis:(M3DVector3f)vector maximum:(NSUInteger*)maximum middle:(NSUInteger*)middle minimum:(NSUInteger*)minimum;

+ (int)nextAxis:(int)axis;

- (void)initializeAccelerometer;

@end
