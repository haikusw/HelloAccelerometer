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
	
	
	
	// :::::::::::::::::::::::::: Box filter :::::::::::::::::::::::::: 
	//	float howmany = ((float)[_samples count]);
	//	m3dLoadVector3f(g, 0.0, 0.0, 0.0);
	//	for (int i = 0; i < gaussianKernalFootprint; i++) {
	//		
	//		g[0] += [[_x_samples objectAtIndex:i] floatValue];
	//		g[1] += [[_y_samples objectAtIndex:i] floatValue];
	//		g[2] += [[_z_samples objectAtIndex:i] floatValue];
	//		
	//	}
	//	m3dScaleVector3f(g, 1.0/howmany);
	
	
	
	
	// :::::::::::::::::::::::::: Gaussian filter :::::::::::::::::::::::::: 	
	m3dLoadVector3f(vec,	0.0, 0.0, 0.0);
	for (int i = 0; i < gaussianKernalFootprint; i++) {
		
		float wgt = gaussianKernalDiscrete(i);
		
		vec[0] += wgt * [[_x_samples objectAtIndex:i] floatValue];
		vec[1] += wgt * [[_y_samples objectAtIndex:i] floatValue];
		vec[2] += wgt * [[_z_samples objectAtIndex:i] floatValue];

	}

	M3DVector3f sph;
	m3dLoadVector3f(sph, 0.0, 0.0, 0.0);
	TIESpherical(vec, sph);
	
	// :::::::::::::::::::::::::: Median filter :::::::::::::::::::::::::: 
	//	m3dLoadVector3f(g, 0.0, 0.0, 0.0);
	//	[_x_samples sortUsingSelector:@selector(compare:)];
	//	[_y_samples sortUsingSelector:@selector(compare:)];
	//	[_z_samples sortUsingSelector:@selector(compare:)];
	//	
	//	int index	= ([_x_samples count] - 1)/2;
	//	g[0] = [[_x_samples objectAtIndex:index] floatValue];
	//	g[1] = [[_y_samples objectAtIndex:index] floatValue];
	//	g[2] = [[_z_samples objectAtIndex:index] floatValue];
	
	
	
	
	
	// Dequeue and discard oldest sample from the queue
	[_x_samples removeObjectAtIndex:0];
	[_y_samples removeObjectAtIndex:0];
	[_z_samples removeObjectAtIndex:0];
	
	// Enqueue most current sample
	[_x_samples addObject:xx];
	[_y_samples addObject:yy];
	[_z_samples addObject:zz];
	
	// We don't care about the magnitude of G, just it's direction
//	m3dNormalizeVectorf(vec);
	
	// Construct transformation matrix
	
	static BOOL gPastSet = NO;
	if (gPastSet == NO) {
		
		m3dCopyVector3f(gPast,			vec);
		m3dCopyVector3f(gSphericalPast,	sph);
		
		gPastSet = YES;
		
		return;
	}
	
	
	float dotProduct	= m3dDotProductf(vec, gPast);	
	float angleBetween	= m3dRadToDeg( acosf( dotProduct / ( m3dGetVectorLengthf(vec) * m3dGetVectorLengthf(gPast) ) ) );
//	NSLog(@"Angle between g(%f %f %f) and gPast(%f %f %f) is %f.", vec[0], vec[1], vec[2], gPast[0], gPast[1], gPast[2], angleBetween);
	
//	float deltaTheta	= fabsf(m3dRadToDeg(sph[1] - gSphericalPast[1]));
//	float deltaPhi		= fabsf(m3dRadToDeg(sph[2] - gSphericalPast[2]));
	
	// Only look at delta phi! Delta theta is highly unstable since it is f(x, y)
	// and x and y are tiny and unstable when the device is flat on a desk.

	static float threshold = 1.0;
	if (angleBetween < threshold) {
		
//		NSLog(@"angleBetween(%f) < threshold(%f) ... BAILING", angleBetween, threshold);
		return;
	} else {
//		NSLog(@"angleBetween(%f) < threshold(%f) ... PROCEEDING", angleBetween, threshold);
	}

	
	m3dCopyVector3f(g,			vec);
	m3dCopyVector3f(gSpherical,	sph);
	
	float dx = fabsf(g[0] - gPast[0]);
	float dy = fabsf(g[1] - gPast[1]);
	
	
	if (dy > dx) {
		
		if (gPast[1] > g[1]) {
			
			m3dCrossProductf(x, g, gPast);	
		} else {
			
			m3dCrossProductf(x, gPast, g);	
		}
		
		m3dCrossProductf(y, x, g);
		
	} else {
		
		if (gPast[0] < g[0]) {
			
			m3dCrossProductf(y, g, gPast);	
		} else {
			
			m3dCrossProductf(y, gPast, g);	
		}
		
		m3dCrossProductf(x, g, y);
		
	}
	
	m3dNormalizeVectorf(x);
	m3dNormalizeVectorf(y);
	
	m3dCrossProductf(z, x, y);	


	
//	// The computed frame is the OpenGL Model View transformation. Store it away
//	m3dLoadIdentity44f(_openGLCameraInverseTransform);
//	MatrixElement(_openGLCameraInverseTransform, 0, 0) = x[0];
//	MatrixElement(_openGLCameraInverseTransform, 1, 0) = x[1];
//	MatrixElement(_openGLCameraInverseTransform, 2, 0) = x[2];
//	
//	MatrixElement(_openGLCameraInverseTransform, 0, 1) = y[0];
//	MatrixElement(_openGLCameraInverseTransform, 1, 1) = y[1];
//	MatrixElement(_openGLCameraInverseTransform, 2, 1) = y[2];
//	
//	MatrixElement(_openGLCameraInverseTransform, 0, 2) = z[0];
//	MatrixElement(_openGLCameraInverseTransform, 1, 2) = z[1];
//	MatrixElement(_openGLCameraInverseTransform, 2, 2) = z[2];
//	
//	// Build the upper 3x3 of the OpenGL style "view" transformation from the transpose of the camera orientation
//	// This is the inversion process. Since these 3x3 matrices are orthonormal a transpose is sufficient to invert
//	m3dLoadIdentity44f(_cameraTransform);	
//	for (int i = 0; i < 3; i++) {
//		for (int j = 0; j < 3; j++) {
//			MatrixElement(_cameraTransform, i, j) = MatrixElement(_openGLCameraInverseTransform, j, i);
//		}
//	}

	
	
	
	
	
	nx_label.text	= [NSString stringWithFormat:@"%.2f", x[0]];
	ny_label.text	= [NSString stringWithFormat:@"%.2f", x[1]];
	nz_label.text	= [NSString stringWithFormat:@"%.2f", x[2]];
	
	ox_label.text	= [NSString stringWithFormat:@"%.2f", y[0]];
	oy_label.text	= [NSString stringWithFormat:@"%.2f", y[1]];
	oz_label.text	= [NSString stringWithFormat:@"%.2f", y[2]];
	
	ax_label.text	= [NSString stringWithFormat:@"%.2f", z[0]];
	ay_label.text	= [NSString stringWithFormat:@"%.2f", z[1]];
	az_label.text	= [NSString stringWithFormat:@"%.2f", z[2]];
	
    rho_label.text			= [NSString stringWithFormat:@"%.2f",	sph[0]				];
    theta_label.text		= [NSString stringWithFormat:@"%.2f",	m3dRadToDeg(sph[1])	];
    phi_label.text			= [NSString stringWithFormat:@"%.2f",	m3dRadToDeg(sph[2])	];
		
	
	

	
	
	
	
	// Update G
	m3dCopyVector3f(gPast,			g);
	m3dCopyVector3f(gSphericalPast,	gSpherical);
	
	
}

#ifdef DEPRICATED
- (void)accelerometer:(UIAccelerometer*)accelerometer didAccelerate:(UIAcceleration*)acceleration {
	
	// Extract G using low pass filter
	m3dCopyVector3f(gPast,          g);
	m3dCopyVector3f(gSphericalPast, gSpherical);
	
	// Raw data
	acc_raw[0] =	acceleration.x;
	acc_raw[1] =	acceleration.y;
	acc_raw[2] =	acceleration.z;

	// unoptimized
	// g = acc_raw * kFilteringFactor + gPast * (1.0 - kFilteringFactor);
	
	// optimized
	// g = kFilteringFactor * (acc_raw - gPast) + gPast;
	
	g[0] = kFilteringFactor * (acc_raw[0] - gPast[0]) + gPast[0];
	g[1] = kFilteringFactor * (acc_raw[1] - gPast[1]) + gPast[1];
	g[2] = kFilteringFactor * (acc_raw[2] - gPast[2]) + gPast[2];
	
	TIESpherical(g, gSpherical);	

	
	// Instantaneous acceleration (G subtracted)
	m3dSubtractVectors3f(acc_instantaneous, acc_raw, g);
	
	M3DVector3f gg_normalized;
	m3dCopyVector3f(gg_normalized,      g);
	m3dNormalizeVectorf(gg_normalized);
	
	// Compute deltas in spherical coordinates
	m3dSubtractVectors3f(deltaSpherical, gSpherical, gSphericalPast);
	
	static float dRho	= 0.0;
	static float dTheta = 0.0;
	static float dPhi	= 0.0;
	
	static int trigger = 3;
	
	if (++trigger > 3) {
		
		dRho	= (  dRho <             fabsf(deltaSpherical[0])  ) ?             fabsf(deltaSpherical[0])  : dRho;
		dTheta = ( dTheta < m3dRadToDeg(fabsf(deltaSpherical[1])) ) ? m3dRadToDeg(fabsf(deltaSpherical[1])) : dTheta;
		dPhi	= (  dPhi < m3dRadToDeg(fabsf(deltaSpherical[2])) ) ? m3dRadToDeg(fabsf(deltaSpherical[2])) : dPhi;
	}
	
	float dx = fabsf(g[0] - gPast[0]);
	float dy = fabsf(g[1] - gPast[1]);

	
	if (dy > dx) {
		
		if (gPast[1] > g[1]) {
			
			m3dCrossProductf(x, g, gPast);	
		} else {
			
			m3dCrossProductf(x, gPast, g);	
		}
		
		m3dCrossProductf(y, x, g);
		
	} else {
		
		if (gPast[0] < g[0]) {
			
			m3dCrossProductf(y, g, gPast);	
		} else {
			
			m3dCrossProductf(y, gPast, g);	
		}
		
		m3dCrossProductf(x, g, y);
		
	}
	
	m3dNormalizeVectorf(x);
	m3dNormalizeVectorf(y);
	
	m3dCrossProductf(z, x, y);	
	
	nx_label.text	= [NSString stringWithFormat:@"%.2f", x[0]];
	ny_label.text	= [NSString stringWithFormat:@"%.2f", x[1]];
	nz_label.text	= [NSString stringWithFormat:@"%.2f", x[2]];
	
	ox_label.text	= [NSString stringWithFormat:@"%.2f", y[0]];
	oy_label.text	= [NSString stringWithFormat:@"%.2f", y[1]];
	oz_label.text	= [NSString stringWithFormat:@"%.2f", y[2]];
	
    ax_label.text	= [NSString stringWithFormat:@"%.2f", z[0]];
    ay_label.text	= [NSString stringWithFormat:@"%.2f", z[1]];
    az_label.text	= [NSString stringWithFormat:@"%.2f", z[2]];
	
    g_x_label.text	= [NSString stringWithFormat:@"%.2f", gg_normalized[0]];
	g_y_label.text	= [NSString stringWithFormat:@"%.2f", gg_normalized[1]];
    g_z_label.text	= [NSString stringWithFormat:@"%.2f", gg_normalized[2]];
	
    delta_rho_label.text	= [NSString stringWithFormat:@"%.2f", dRho];
	delta_theta_label.text	= [NSString stringWithFormat:@"%.2f", dTheta];
    delta_phi_label.text	= [NSString stringWithFormat:@"%.2f", dPhi];
	
}
#endif

@end
