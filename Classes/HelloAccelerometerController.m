//
//  HelloAccelerometerController.h
//  HelloAccelerometer
//
//  Created by turner on 5/6/09.
//  Copyright Douglass Turner Consulting 2009. All rights reserved.
//


#import "AccelerometerSimulation.h"
#import "HelloAccelerometerController.h"

#define kAccelerometerFrequency		(30)
#define kFilteringFactor			(0.1)

static BOOL kernalLookuptableHasBeenBuilt = NO;
NSMutableArray* kernalLookupTable;
float gaussianKernalHelper(float in);
float gaussianKernal(float in);
float gaussianKernalDiscrete(int index);

@implementation HelloAccelerometerController

- (void)initializeAccelerometer {
	
	_x_samples		= [[NSMutableArray alloc] init];
	_y_samples		= [[NSMutableArray alloc] init];
	_z_samples		= [[NSMutableArray alloc] init];

	kernalLookupTable = [[NSMutableArray alloc] init];
	
	[[UIAccelerometer sharedAccelerometer] setUpdateInterval:(1.0 / kAccelerometerFrequency)];
	[[UIAccelerometer sharedAccelerometer] setDelegate:self];
	
	// !!!!!!!!!!!!!!! DISABLE IDLE SCREEN !!!!!!!!!!!!!!!!!!
	[UIApplication sharedApplication].idleTimerDisabled = YES;
	
}

const int gaussianKernalFootprint		= 15;
const float gaussianKernalScaleFactor	= 2.70;
const float gaussianKernalSigma			= 1.0;

static float gaussianKernalNormalization	= 0.0;

float gaussianKernalHelper(float in) {
	
	// in is assumed to be [-gaussianKernalScaleFactor < in < gaussianKernalScaleFactor]
	float g = expf(-(in * in)/(2.0 * gaussianKernalSigma * gaussianKernalSigma)) / (sqrtf(2.0 * M_PI) * gaussianKernalSigma);
	
	return g;
}

// index is assumed to be [0 <= index <= gaussianKernalFootprint - 1]
float gaussianKernalDiscrete(int index) {
	
	if (kernalLookuptableHasBeenBuilt == YES) {
		return [[kernalLookupTable objectAtIndex:index] floatValue] / gaussianKernalNormalization;
	}
	
	int delta = ((gaussianKernalFootprint - 1)/2);
	float w	= ((float)(index - delta))/((float)delta);
	
	return gaussianKernal( w );

}

// in is assumed to be [-1 <= in <= 1]
float gaussianKernal(float in) {

	static BOOL gaussianKernalHasBeenNormalized	= NO;
	if (gaussianKernalHasBeenNormalized == NO) {
		
		int i;
		float f;
		
		int howmany		= (gaussianKernalFootprint - 1) / 2;
		float increment	= gaussianKernalScaleFactor / ((float) howmany);
		
		for (i = 0, f = -(increment * ((float)howmany)); i < gaussianKernalFootprint; i++, f += increment) {
			
			// gaussianKernalHelper(g) takes [-gaussianKernalScaleFactor < g < gaussianKernalScaleFactor]
			gaussianKernalNormalization += gaussianKernalHelper(f);
		}
						
		// Sanity check that total weight under kernal is 1
//		float sum = 0.0;
//		for (i = 0, f = -(increment * ((float)howmany)); i < gaussianKernalFootprint; i++, f += increment) {
//			
//			// gaussianKernalHelper(g) takes [-gaussianKernalScaleFactor < g < gaussianKernalScaleFactor]
//			sum += gaussianKernalHelper(f);
//		}
		
		gaussianKernalHasBeenNormalized	= YES;
	}

	if (kernalLookuptableHasBeenBuilt == NO) {
		
		int i;
		float f;
		
		int howmany		= (gaussianKernalFootprint - 1) / 2;
		float increment	= gaussianKernalScaleFactor / ((float) howmany);
		
		// Build lookup table
		for (i = 0, f = -(increment * ((float)howmany)); i < gaussianKernalFootprint; i++, f += increment) {
			
			float wgt = gaussianKernalHelper(f);
			[kernalLookupTable addObject:[NSNumber numberWithFloat:wgt]];
		}
		
		kernalLookuptableHasBeenBuilt = YES;
	}
	
	float g = gaussianKernalHelper(in * gaussianKernalScaleFactor);
	
	return g / gaussianKernalNormalization;
}

// Build x and y axes. We use a simple heuristic: observe the deltas in x and y,
// which ever is greater indicates the axis to create. Greater delta in x implies
// calculation of the y-axis by z X x. Greater delta in y implies calculation of
// the x-axis by y X z. Simple, but it appears to work - 9 May 2009
//
//
// The resultant coordinate frame is orientated as follows
//
// x-axis to the right
//
// y-axis upwards
//
// z-axis normal to the device and facing outwards
//        towards the user
//
//      +---------------+
//		| +-----------+ |
//		| |		Y     | |
//		| |		|     | |
//		| |		|     | |
//		| |		|     | |
//		| |		+----------X
//		| |		      | |
//		| |			  | |
//		| |		      | |
//		| |		      | |
//		| +-----------+ |
//		|	  +---+	    |
//		|     | O |	    |
//		|	  +---+	    |
//      +---------------+
//
// UIAccelerometerDelegate method
- (void)accelerometer:(UIAccelerometer*)accelerometer didAccelerate:(UIAcceleration*)acceleration {
	
	// most current raw sample
	NSNumber *xx = [NSNumber numberWithFloat:acceleration.x];
	NSNumber *yy = [NSNumber numberWithFloat:acceleration.y];
	NSNumber *zz = [NSNumber numberWithFloat:acceleration.z];
	
	M3DVector3f vec;
	m3dLoadVector3f(vec, [xx floatValue], [yy floatValue], [zz floatValue]);
	
	// Accumulate raw samples
	if ([_x_samples count] < gaussianKernalFootprint) {
		
		// Enqueue
		[_x_samples addObject:xx];
		[_y_samples addObject:yy];
		[_z_samples addObject:zz];
		
		return;
	}
	
	// :::::::::::::::::::::::::: Gaussian filter :::::::::::::::::::::::::: 	
	m3dLoadVector3f(vec,	0.0, 0.0, 0.0);
	for (int i = 0; i < gaussianKernalFootprint; i++) {
		
		float wgt = gaussianKernalDiscrete(i);
		
		vec[0] += wgt * [[_x_samples objectAtIndex:i] floatValue];
		vec[1] += wgt * [[_y_samples objectAtIndex:i] floatValue];
		vec[2] += wgt * [[_z_samples objectAtIndex:i] floatValue];

	}

	
	// Dequeue and discard oldest sample from the queue
	[_x_samples removeObjectAtIndex:0];
	[_y_samples removeObjectAtIndex:0];
	[_z_samples removeObjectAtIndex:0];
	
	// Enqueue most current sample
	[_x_samples addObject:xx];
	[_y_samples addObject:yy];
	[_z_samples addObject:zz];
	
	// Do nothing until we have a past value of "G"
	static BOOL gPastSet = NO;
	if (gPastSet == NO) {
		
		m3dCopyVector3f(gPast, vec);		
		gPastSet = YES;
		
		return;
	}
	
	
	// Establish a threshold below which we discard samples. We need a large enough delta (ie, dot product) to compute a valid
	// coordinate frame
	float dotProduct	= m3dDotProductf(vec, gPast);	
	float angleBetween	= m3dRadToDeg( acosf( dotProduct / ( m3dGetVectorLengthf(vec) * m3dGetVectorLengthf(gPast) ) ) );
//	NSLog(@"Angle between g(%f %f %f) and gPast(%f %f %f) is %f.", vec[0], vec[1], vec[2], gPast[0], gPast[1], gPast[2], angleBetween);
		
	static float threshold = 1.0;
	if (angleBetween < threshold) {
		
//		NSLog(@"angleBetween(%f) < threshold(%f) ... BAILING", angleBetween, threshold);
		return;
	} else {
//		NSLog(@"angleBetween(%f) < threshold(%f) ... PROCEEDING", angleBetween, threshold);
	}
	
	m3dCopyVector3f(g, vec);
	TIESpherical(g, gSpherical);
	
	M3DVector3f g_yzPlane, gPast_yzPlane;
	m3dLoadVector3f(    g_yzPlane,	0.0,     g[1],     g[2]);
	m3dLoadVector3f(gPast_yzPlane,	0.0, gPast[1], gPast[2]);

//	float yz_dotProduct			= m3dDotProductf(g_yzPlane, gPast_yzPlane);	
//	float yz_angle_between	= m3dRadToDeg( acosf( yz_dotProduct / ( m3dGetVectorLengthf(g_yzPlane) * m3dGetVectorLengthf(gPast_yzPlane) ) ) );

	M3DVector3f g_xzPlane, gPast_xzPlane;
	m3dLoadVector3f(    g_xzPlane,	    g[0], 0.0,     g[2]);
	m3dLoadVector3f(gPast_xzPlane,	gPast[0], 0.0, gPast[2]);
	
//	float xz_dotProduct			= m3dDotProductf(g_xzPlane, gPast_xzPlane);	
//	float xz_angle_between	= m3dRadToDeg( acosf( xz_dotProduct / ( m3dGetVectorLengthf(g_xzPlane) * m3dGetVectorLengthf(gPast_xzPlane) ) ) );
	
//	if (yz_dotProduct > xz_dotProduct) {
//		NSLog(@"YZ: DotProduct(%.3f) > XZ: DotProduct(%.3f)", yz_dotProduct, xz_dotProduct);
//	} else {
//		NSLog(@"YZ: DotProduct(%.3f) < XZ: DotProduct(%.3f)", yz_dotProduct, xz_dotProduct);
//	}

	float dx = fabsf(g[0] - gPast[0]);
	float dy = fabsf(g[1] - gPast[1]);
	
	if (dy > dx) {
//	if (yz_dotProduct > xz_dotProduct) {
		
		if (gPast[1] > g[1]) {
			
			m3dCrossProductf(exe, g, gPast);	
		} else {
			
			m3dCrossProductf(exe, gPast, g);	
		}
		
		m3dCrossProductf(wye, exe, g);
		
	} else {
		
		if (gPast[0] < g[0]) {
			
			m3dCrossProductf(wye, g, gPast);	
		} else {
			
			m3dCrossProductf(wye, gPast, g);	
		}
		
		m3dCrossProductf(exe, g, wye);
		
	}
	
	m3dNormalizeVectorf(exe);
	m3dNormalizeVectorf(wye);
	
	m3dCrossProductf(zee, exe, wye);	


	
	// The computed frame is the OpenGL Model View transformation. Store it away
	m3dLoadIdentity44f(openGLModelViewTransform);
	MatrixElement(openGLModelViewTransform, 0, 0) = exe[0];
	MatrixElement(openGLModelViewTransform, 1, 0) = exe[1];
	MatrixElement(openGLModelViewTransform, 2, 0) = exe[2];
	
	MatrixElement(openGLModelViewTransform, 0, 1) = wye[0];
	MatrixElement(openGLModelViewTransform, 1, 1) = wye[1];
	MatrixElement(openGLModelViewTransform, 2, 1) = wye[2];
	
	MatrixElement(openGLModelViewTransform, 0, 2) = zee[0];
	MatrixElement(openGLModelViewTransform, 1, 2) = zee[1];
	MatrixElement(openGLModelViewTransform, 2, 2) = zee[2];
	
	// Build the upper 3x3 of the OpenGL style "view" transformation from the transpose of the camera orientation
	// This is the inversion process. Since these 3x3 matrices are orthonormal a transpose is sufficient to invert
	m3dLoadIdentity44f(modelingTransform);	
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			MatrixElement(modelingTransform, i, j) = MatrixElement(openGLModelViewTransform, j, i);
		}
	}

	
	
	
	
	
	nx_label.text	= [NSString stringWithFormat:@"%.2f", exe[0]];
	ny_label.text	= [NSString stringWithFormat:@"%.2f", exe[1]];
	nz_label.text	= [NSString stringWithFormat:@"%.2f", exe[2]];
	
	ox_label.text	= [NSString stringWithFormat:@"%.2f", wye[0]];
	oy_label.text	= [NSString stringWithFormat:@"%.2f", wye[1]];
	oz_label.text	= [NSString stringWithFormat:@"%.2f", wye[2]];
	
	ax_label.text	= [NSString stringWithFormat:@"%.2f", zee[0]];
	ay_label.text	= [NSString stringWithFormat:@"%.2f", zee[1]];
	az_label.text	= [NSString stringWithFormat:@"%.2f", zee[2]];
	
      rho_label.text	= [NSString stringWithFormat:@"%.2f",	gSpherical[0]				];
    theta_label.text	= [NSString stringWithFormat:@"%.1f",	m3dRadToDeg(gSpherical[1])	];
      phi_label.text	= [NSString stringWithFormat:@"%.1f",	m3dRadToDeg(gSpherical[2])	];
		
	
	

	
	
	
	
	// Update G
	m3dCopyVector3f(gPast, g);

	
}

@end
