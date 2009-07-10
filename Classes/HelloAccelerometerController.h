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
	M3DVector3f				exe;
	
    IBOutlet UILabel		*ox_label;
    IBOutlet UILabel		*oy_label;
    IBOutlet UILabel		*oz_label;
	M3DVector3f				wye;
	
    IBOutlet UILabel		*ax_label;
    IBOutlet UILabel		*ay_label;
    IBOutlet UILabel		*az_label;
	M3DVector3f				zee;
	
	M3DMatrix44f			modelingTransform;
	M3DMatrix44f			openGLModelViewTransform;
	
	// Spherical
    IBOutlet UILabel		*rho_label;
    IBOutlet UILabel		*theta_label;
    IBOutlet UILabel		*phi_label;
	M3DVector3f				gSpherical;
	
	// Instantaneous acceleration
    IBOutlet UILabel		*acc_instantaneous_x_label;
    IBOutlet UILabel		*acc_instantaneous_y_label;
    IBOutlet UILabel		*acc_instantaneous_z_label;
	M3DVector3f				acc_instantaneous;
	
	NSMutableArray* _x_samples;
	NSMutableArray* _y_samples;
	NSMutableArray* _z_samples;
}

+(NSString*)signPrint:(float)f;

+ (NSString*)signedAxisName:(M3DVector3f)vector axis:(NSUInteger)axis;

+(float)sign:(float)f;

+ (void)dominantAxis:(M3DVector3f)vector maximum:(NSUInteger*)maximum middle:(NSUInteger*)middle minimum:(NSUInteger*)minimum;

+ (int)nextAxis:(int)axis;

- (void)initializeAccelerometer;

@end
