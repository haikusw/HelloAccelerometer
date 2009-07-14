//
//  HelloAccelerometerController.h
//  HelloAccelerometer
//
//  Created by turner on 5/6/09.
//  Copyright Douglass Turner Consulting 2009. All rights reserved.
//


#import "AccelerometerSimulation.h"
#import "HelloAccelerometerController.h"


static M3DVector3f _compensationVector = { -0.06, 0.14, 0.08 };
//static M3DVector3f _xAxis = { 1.0, 0.0, 0.0 };
//static M3DVector3f _yAxis = { 0.0, 1.0, 0.0 };
//static M3DVector3f _zAxis = { 0.0, 0.0, 1.0 };

#define kAccelerometerFrequency			(60)
#define kFilteringFactor				(0.1)
#define kIgnoreGravitySampleThreshold	(360.0/200.0)

static BOOL kernalLookuptableHasBeenBuilt = NO;
NSMutableArray* kernalLookupTable;
float gaussianKernalHelper(float in);
float gaussianKernal(float in);
float gaussianKernalDiscrete(int index);

@implementation HelloAccelerometerController

+(NSString*)signPrint:(float)f {
	
	if (f < 0.0) return @"";
	
	return @"+";
}

+ (NSString*)signedAxisName:(M3DVector3f)vector axis:(NSUInteger)axis {
	
	NSString* result = nil;
	
	switch (axis) {
			
		case 0:
			result = @"+X";
			if (vector[axis] < 0.0) {
				result = @"-X";
			}
			break;
		case 1:
			result = @"+Y";
			if (vector[axis] < 0.0) {
				result = @"-Y";
			}
			break;
		case 2:
			result = @"+Z";
			if (vector[axis] < 0.0) {
				result = @"-Z";
			}
			break;
		default:
			result = @"Huh?";
	}
	
	return result;
};

+(float)sign:(float)f {
	
	if (f < 0.0f) return -1.0f;
	
	return 1.0f;
}

+ (void)dominantAxis:(M3DVector3f)vector maximum:(NSUInteger*)maximum middle:(NSUInteger*)middle minimum:(NSUInteger*)minimum {

	M3DVector3f fabulous;
	m3dLoadVector3f(fabulous, fabsf(vector[0]), fabsf(vector[1]), fabsf(vector[2]));
	
	
	if (fabulous[0] > fabulous[1] && fabulous[0] > fabulous[2]) {
		
		*maximum = 0;
		
		if (fabulous[1] > fabulous[2]) {
			*middle		= 1;
			*minimum	= 2;
		} else {
			*middle		= 2;
			*minimum	= 1;
		}
		
	}
	
	if (fabulous[1] > fabulous[2] && fabulous[1] > fabulous[0]) {
		
		*maximum = 1;
		
		if (fabulous[2] > fabulous[0]) {
			*middle		= 2;
			*minimum	= 0;
		} else {
			*middle		= 0;
			*minimum	= 2;
		}
		
	}
	
	if (fabulous[2] > fabulous[0] && fabulous[2] > fabulous[1]) {
		
		*maximum = 2;
		
		if (fabulous[0] > fabulous[1]) {
			*middle		= 0;
			*minimum	= 1;
		} else {
			*middle		= 1;
			*minimum	= 0;
		}
		
	}
	
}

+ (int)nextAxis:(int)axis {
	
	int next;
	
	next = (axis + 1) % 3;
	
	return next;
}

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

static const int gaussianKernalFootprint		= 15;
static const float gaussianKernalScaleFactor	= 2.70;
static const float gaussianKernalSigma			= 1.0;

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

- (void)accelerometer:(UIAccelerometer*)accelerometer didAccelerate:(UIAcceleration*)acceleration {

	
	
	// most current raw sample
	NSNumber *xx = [NSNumber numberWithFloat:acceleration.x];
	NSNumber *yy = [NSNumber numberWithFloat:acceleration.y];
	NSNumber *zz = [NSNumber numberWithFloat:acceleration.z];

	
	
	
	// Accumulate a sufficient numbers samples before going on to calculate "G".
	if ([_x_samples count] < gaussianKernalFootprint) {
		
		// Enqueue
		[_x_samples addObject:xx];
		[_y_samples addObject:yy];
		[_z_samples addObject:zz];
		
		return;
	}

	
	
	
	// Using a Gaussian kernal, create a single filtered sample from the raw samples.
	M3DVector3f vec;
	m3dLoadVector3f(vec,	0.0, 0.0, 0.0);
	for (int i = 0; i < gaussianKernalFootprint; i++) {
		
		float wgt = gaussianKernalDiscrete(i);
		
		vec[0] += wgt * [[_x_samples objectAtIndex:i] floatValue];
		vec[1] += wgt * [[_y_samples objectAtIndex:i] floatValue];
		vec[2] += wgt * [[_z_samples objectAtIndex:i] floatValue];

	}


	
	
	// Dequeue and discard the oldest sample
	[_x_samples removeObjectAtIndex:0];
	[_y_samples removeObjectAtIndex:0];
	[_z_samples removeObjectAtIndex:0];

	
	
	
	
	// Enqueue the most current sample
	[_x_samples addObject:xx];
	[_y_samples addObject:yy];
	[_z_samples addObject:zz];
	


	
	
	// Add calibration compensation
	m3dAddVectors2f(vec, vec, _compensationVector);

	
	
	
	// Do nothing until we have a past value of "G"
	static BOOL isGPastSet = NO;
	if (isGPastSet == NO) {
		
		// The most recent sample becomes gPast
		m3dCopyVector3f(gPast, vec);
		isGPastSet = YES;
		
		return;
	}

	
	// If the angle between "g" and "gPast" is below the threshold, bail. We need a "g" and "gPast" that are sufficiently
	// separated in space to accurately calculate a new coordinate frame. We will keep "gPast" as we examine candidates
	// for "g".
	float dotProduct	= m3dDotProductf(vec, gPast);	
	float angleBetween	= m3dRadToDeg( acosf( dotProduct / ( m3dGetVectorLengthf(vec) * m3dGetVectorLengthf(gPast) ) ) );

	if (angleBetween < kIgnoreGravitySampleThreshold) {
		
		return;
	}


	
	
	// Accept vec as current value of g
	m3dCopyVector3f(g, vec);
	
	
	// Determine dominant axis pattern for g and gPast
	NSUInteger a, b, c;
	[HelloAccelerometerController dominantAxis:g maximum:&a middle:&b minimum:&c];
	
	NSUInteger aPast, bPast, cPast;
	[HelloAccelerometerController dominantAxis:gPast maximum:&aPast middle:&bPast minimum:&cPast];
	
	// We want the axis dominance pattern for "g" and "gPast" to be identical. We don't care about
	// the magnitude or sign of corresponding axes, just that they follow the same pattern of precedence.
	if (a != aPast) {
		goto updateGPast;
	}
	
	if (b != bPast) {
		goto updateGPast;
	}
	
	if (c != cPast) {
		goto updateGPast;
	}

	
	// The component with least displacement corresponds to the axis (X, Y, or Z) we will calculate from
	// the cross product of "g" and "gPast"
	M3DVector3f deltaVector;
	m3dSubtractVectors3f(deltaVector, g, gPast);
	
	NSUInteger aDelta, bDelta, cDelta;
	[HelloAccelerometerController dominantAxis:deltaVector maximum:&aDelta middle:&bDelta minimum:&cDelta];

	
	
	M3DVector3f deltaMult;
	m3dCopyVector3f(deltaMult, g);
	deltaMult[0] *= deltaVector[0];
	deltaMult[1] *= deltaVector[1];
	deltaMult[2] *= deltaVector[2];
	
	NSUInteger aDeltaMult, bDeltaMult, cDeltaMult;
	[HelloAccelerometerController dominantAxis:deltaMult maximum:&aDeltaMult middle:&bDeltaMult minimum:&cDeltaMult];

	
	
	
	// Bail if the dominant delta is along the z-component
	if (cDeltaMult == 2) {
//		NSLog(@"WARNING! WARNING! Z-COMPONENT OF G HAS MINIMUM DELTA. BAILING ...");
		goto updateGPast;
	} // if (cDeltaMult == 2) 
	
		

	
	TIESpherical(g, gSpherical);
	
	// Offset phi to a more convenient value.
	gSpherical[2] = M_PI - gSpherical[2];

	
	
	//	float dot_x_axis				= m3dDotProductf(g, _xAxis);  
	//	float dot_x_axis_angle_between	= m3dRadToDeg( acosf( dot_x_axis / ( m3dGetVectorLengthf(g) * m3dGetVectorLengthf(_xAxis) ) ) );
	//
	//	float dot_y_axis				= m3dDotProductf(g, _yAxis);  
	//	float dot_y_axis_angle_between	= m3dRadToDeg( acosf( dot_y_axis / ( m3dGetVectorLengthf(g) * m3dGetVectorLengthf(_yAxis) ) ) );
	
	
	
	
	
	M3DVector3f axis;		
	NSString* axisName = @"nuthin'";
	if (cDelta == 0) axisName = @"EXE";
	if (cDelta == 1) axisName = @"WYE";
	if (cDelta == 2) axisName = @"ZEE";
			
	m3dCrossProductf(axis, g,		gPast);				
	m3dNormalizeVectorf(axis);
	
	if ([HelloAccelerometerController sign:axis[ cDeltaMult ]] < 0.0) {
		m3dScaleVector3f(axis, -1.0);
	}
	
	NSString* xs = @"";
	NSString* ys = @"";

	// We have calculated exe (cDelta == 0) from the delta between g and gPast.
	if (cDeltaMult == 0) {
		
		xs = @"*";
		m3dCopyVector3f(exe, axis);
		
		m3dCrossProductf(wye, exe, g);			
		m3dNormalizeVectorf(wye);
		
	} // if (cDeltaMult == 0)
	
	
	
	
	// We have calculated wye (cDelta == 1) from the delta between g and gPast.
	if (cDeltaMult == 1) {
		
		ys = @"*";
		m3dCopyVector3f(wye, axis);
		
		m3dCrossProductf(exe, g, wye);			
		m3dNormalizeVectorf(exe);
				
	} // if (cDeltaMult == 1)
	
	// Finally, zee = exe X wye;
	m3dCrossProductf(zee, exe, wye);			

	
	
	
//	NSLog(@"%@%@: %@%.3f %@%.3f %@%.3f		%@%@: %@%.3f %@%.3f %@%.3f		%@: %@%.3f %@%.3f %@%.3f		G:%@ %@ %@		Delta: %@ %@ %@", 
//		  
//		  xs,
//		  @"EXE", 
//		  [HelloAccelerometerController signPrint:exe[0]],
//		  exe[0], 
//		  [HelloAccelerometerController signPrint:exe[1]],
//		  exe[1], 
//		  [HelloAccelerometerController signPrint:exe[2]],
//		  exe[2], 
//		  
//		  ys,
//		  @"WYE", 
//		  [HelloAccelerometerController signPrint:wye[0]],
//		  wye[0], 
//		  [HelloAccelerometerController signPrint:wye[1]],
//		  wye[1], 
//		  [HelloAccelerometerController signPrint:wye[2]],
//		  wye[2], 
//		  
//		  @"ZEE", 
//		  [HelloAccelerometerController signPrint:zee[0]],
//		  zee[0], 
//		  [HelloAccelerometerController signPrint:zee[1]],
//		  zee[1], 
//		  [HelloAccelerometerController signPrint:zee[2]],
//		  zee[2], 
//		  
//		  [HelloAccelerometerController signedAxisName:g axis:a],
//		  [HelloAccelerometerController signedAxisName:g axis:b],
//		  [HelloAccelerometerController signedAxisName:g axis:c],
//		  
//		  [HelloAccelerometerController signedAxisName:deltaVector axis:aDelta],
//		  [HelloAccelerometerController signedAxisName:deltaVector axis:bDelta],
//		  [HelloAccelerometerController signedAxisName:deltaVector axis:cDelta]
//		  
//		  );
	
	
		
		
	
	
	
	
	// exe, wye, and zee are the basis vectors of the camera transform.
	m3dLoadIdentity44f(cameraTransform);
	MatrixElement(cameraTransform, 0, 0) = exe[0];
	MatrixElement(cameraTransform, 1, 0) = exe[1];
	MatrixElement(cameraTransform, 2, 0) = exe[2];
	
	MatrixElement(cameraTransform, 0, 1) = wye[0];
	MatrixElement(cameraTransform, 1, 1) = wye[1];
	MatrixElement(cameraTransform, 2, 1) = wye[2];
	
	MatrixElement(cameraTransform, 0, 2) = zee[0];
	MatrixElement(cameraTransform, 1, 2) = zee[1];
	MatrixElement(cameraTransform, 2, 2) = zee[2];

	// Inversion of the camera transform yields the modeling transform. Since the camera transform
	// is orthonormal, its inverse is the same as its transpose.
	M3DMatrix44f pureOrientation;
	m3dLoadIdentity44f(pureOrientation);	
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			MatrixElement(pureOrientation, i, j) = MatrixElement(cameraTransform, j, i);
		}
	}

	M3DMatrix44f translation;
	TIEMatrix4x4LoadTranslation(translation, 0.0, 0.0, 16.0);
	
	// TIEMatrix4x4Concatenation(B, A, B*A) assuming pt' results from B * A * pt
	TIEMatrix4x4Concatenation(pureOrientation, translation, modelingTransform);

	
	
	
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!! TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// Concatenate the modeling orientation transform (from above) with a translation vector
	// then invert into the OpenGL camera transform using my code from Teapot Toy.
	//
	// This is the relevant snippet of code. "p" is the camera translation vector
	//	MatrixElement(_openGLCameraInverseTransform, 0, 3) = -m3dDotProductf(p, n);
	//	MatrixElement(_openGLCameraInverseTransform, 1, 3) = -m3dDotProductf(p, o);
	//	MatrixElement(_openGLCameraInverseTransform, 2, 3) = -m3dDotProductf(p, a);

	
	
	nx_label.text	= [NSString stringWithFormat:@"%.2f", MatrixElement(cameraTransform, 0, 0)];
	ny_label.text	= [NSString stringWithFormat:@"%.2f", MatrixElement(cameraTransform, 1, 0)];
	nz_label.text	= [NSString stringWithFormat:@"%.2f", MatrixElement(cameraTransform, 2, 0)];
	
	ox_label.text	= [NSString stringWithFormat:@"%.2f", MatrixElement(cameraTransform, 0, 1)];
	oy_label.text	= [NSString stringWithFormat:@"%.2f", MatrixElement(cameraTransform, 1, 1)];
	oz_label.text	= [NSString stringWithFormat:@"%.2f", MatrixElement(cameraTransform, 2, 1)];
	
	ax_label.text	= [NSString stringWithFormat:@"%.2f", MatrixElement(cameraTransform, 0, 2)];
	ay_label.text	= [NSString stringWithFormat:@"%.2f", MatrixElement(cameraTransform, 1, 2)];
	az_label.text	= [NSString stringWithFormat:@"%.2f", MatrixElement(cameraTransform, 2, 2)];
	
	

	
	nx_model_label.text	= [NSString stringWithFormat:@"%.2f", MatrixElement(modelingTransform, 0, 0)];
	ny_model_label.text	= [NSString stringWithFormat:@"%.2f", MatrixElement(modelingTransform, 1, 0)];
	nz_model_label.text	= [NSString stringWithFormat:@"%.2f", MatrixElement(modelingTransform, 2, 0)];
	
	ox_model_label.text	= [NSString stringWithFormat:@"%.2f", MatrixElement(modelingTransform, 0, 1)];
	oy_model_label.text	= [NSString stringWithFormat:@"%.2f", MatrixElement(modelingTransform, 1, 1)];
	oz_model_label.text	= [NSString stringWithFormat:@"%.2f", MatrixElement(modelingTransform, 2, 1)];
	
	ax_model_label.text	= [NSString stringWithFormat:@"%.2f", MatrixElement(modelingTransform, 0, 2)];
	ay_model_label.text	= [NSString stringWithFormat:@"%.2f", MatrixElement(modelingTransform, 1, 2)];
	az_model_label.text	= [NSString stringWithFormat:@"%.2f", MatrixElement(modelingTransform, 2, 2)];
	
	px_model_label.text	= [NSString stringWithFormat:@"%.2f", MatrixElement(modelingTransform, 0, 3)];
	py_model_label.text	= [NSString stringWithFormat:@"%.2f", MatrixElement(modelingTransform, 1, 3)];
	pz_model_label.text	= [NSString stringWithFormat:@"%.2f", MatrixElement(modelingTransform, 2, 3)];
	
	
	
	
      rho_label.text	= [NSString stringWithFormat:@"%.2f",	gSpherical[0]				];
    theta_label.text	= [NSString stringWithFormat:@"%.1f",	m3dRadToDeg(gSpherical[1])	];
      phi_label.text	= [NSString stringWithFormat:@"%.1f",	m3dRadToDeg(gSpherical[2])	];

	
	
	g_x_label.text	= [NSString stringWithFormat:@"%.2f", g[0]];
	g_y_label.text	= [NSString stringWithFormat:@"%.2f", g[1]];
	g_z_label.text	= [NSString stringWithFormat:@"%.2f", g[2]];
	
	
	

	
	
	
	
	// Update G
updateGPast:
	m3dCopyVector3f(gPast, g);

	
}

@end
