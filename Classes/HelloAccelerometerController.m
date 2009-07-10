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

#define kAccelerometerFrequency		(30)
#define kFilteringFactor			(0.1)

static BOOL kernalLookuptableHasBeenBuilt = NO;
NSMutableArray* kernalLookupTable;
float gaussianKernalHelper(float in);
float gaussianKernal(float in);
float gaussianKernalDiscrete(int index);

@implementation HelloAccelerometerController

+(NSString*)signPrint:(float)f {
	
	if (f < 0.0) return @"-";
	
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
	

	
	// Add calibration compensation
	m3dAddVectors2f(vec, vec, _compensationVector);

	
	
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
	
	static float threshold = 1.0/0.5;
	if (angleBetween < threshold) {
		
		return;
	}

	
	
	
	// Accept vec as current value of g
	m3dCopyVector3f(g, vec);
	
	
	
	// Determine dominant axis pattern for g and gPast
	NSUInteger a, b, c;
	[HelloAccelerometerController dominantAxis:g maximum:&a middle:&b minimum:&c];
	
	NSUInteger aPast, bPast, cPast;
	[HelloAccelerometerController dominantAxis:gPast maximum:&aPast middle:&bPast minimum:&cPast];
	
	// Determine which axis has the least displacement. That is the axis we will
	// calculate.
	M3DVector3f deltaVector;
	m3dSubtractVectors3f(deltaVector, g, gPast);
	
	NSUInteger aDelta, bDelta, cDelta;
	[HelloAccelerometerController dominantAxis:deltaVector maximum:&aDelta middle:&bDelta minimum:&cDelta];
	
	if (a == aPast && b == bPast && c == cPast && cDelta != 2) {

		
		TIESpherical(g, gSpherical);
		
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
		
		if ([HelloAccelerometerController sign:axis[cDelta]] < 0.0) {
			m3dScaleVector3f(axis, -1.0);
		}

		// For now, only handle the case where the device is facing up.
		if (g[a] < 0.0) {
			
			M3DVector3f resultOfCross;
			
			// Handle the case where have calculated the x-axis (cDelta == 0) and the z-component of the G vector is dominant axis (a == 2).
			if (cDelta == 0 && a == 2) {
				
				m3dCrossProductf(resultOfCross, axis, g);			
				m3dNormalizeVectorf(resultOfCross);
				
				NSLog(@"%@: %.3f %.3f %.3f	CrossingAxis: %.3f %.3f %.3f	ResultOfCross: %.3f %.3f %.3f", 
					  axisName, 
					  axis[0], axis[1], axis[2], 
					  g[0], g[1], g[2],
					  resultOfCross[0], resultOfCross[1], resultOfCross[2]
					  );
				
			} // if (cDelta == 0 && a == 2)
			
			// Handle the case where have calculated the y-axis (cDelta == 1) and the z-component of the G vector is dominant axis (a == 2).
			if (cDelta == 1 && a == 2) {
				
				m3dCrossProductf(resultOfCross, g, axis);			
				m3dNormalizeVectorf(resultOfCross);
				
				NSLog(@"%@: %.3f %.3f %.3f	CrossingAxis: %.3f %.3f %.3f	ResultOfCross: %.3f %.3f %.3f", 
					  axisName, 
					  axis[0], axis[1], axis[2], 
					  g[0], g[1], g[2],
					  resultOfCross[0], resultOfCross[1], resultOfCross[2]
					  );
				
			} // if (cDelta == 1 && a == 2)
			
		}
	
		
		
		
		
//		if (cDelta == 1 && a == 2) {
//			m3dCrossProductf(resultOfCross, axis, g);			
//		}
		
//		NSLog(@"%@: %.3f %.3f %.3f		G: %@ %.3f %@ %.3f %@ %.3f		GPAST: %@ %.3f %@ %.3f %@ %.3f		Delta: %@ %.3f %@ %.3f %@ %.3f",
//			  axisName,
//			  axis[0], 
//			  axis[1], 
//			  axis[2], 
//			  
//			  [HelloAccelerometerController signedAxisName:g axis:a], 
//			  g[a],    
//			  [HelloAccelerometerController signedAxisName:g axis:b], 
//			  g[b],    
//			  [HelloAccelerometerController signedAxisName:g axis:c],
//			  g[c],
//			  
//			  [HelloAccelerometerController signedAxisName:gPast axis:aPast], 
//			  gPast[aPast],    
//			  [HelloAccelerometerController signedAxisName:gPast axis:bPast], 
//			  gPast[bPast],    
//			  [HelloAccelerometerController signedAxisName:gPast axis:cPast],
//			  gPast[cPast],
//			  
//			  [HelloAccelerometerController signedAxisName:deltaVector axis:aDelta],
//			  deltaVector[aDelta], 
//			  [HelloAccelerometerController signedAxisName:deltaVector axis:bDelta],
//			  deltaVector[bDelta], 
//			  [HelloAccelerometerController signedAxisName:deltaVector axis:cDelta],
//			  deltaVector[cDelta]
//			  );

		
		
		
		
		
		
		
		
		
		
		
		
	} // if (a == aPast && b == bPast && c == cPast && cDelta != 2)
	
	
	
	
	
	// The computed frame is the OpenGL Model View transformation. Store it away
//	m3dLoadIdentity44f(openGLModelViewTransform);
//	MatrixElement(openGLModelViewTransform, 0, 0) = exe[0];
//	MatrixElement(openGLModelViewTransform, 1, 0) = exe[1];
//	MatrixElement(openGLModelViewTransform, 2, 0) = exe[2];
//	
//	MatrixElement(openGLModelViewTransform, 0, 1) = wye[0];
//	MatrixElement(openGLModelViewTransform, 1, 1) = wye[1];
//	MatrixElement(openGLModelViewTransform, 2, 1) = wye[2];
//	
//	MatrixElement(openGLModelViewTransform, 0, 2) = zee[0];
//	MatrixElement(openGLModelViewTransform, 1, 2) = zee[1];
//	MatrixElement(openGLModelViewTransform, 2, 2) = zee[2];
	
	// Build the upper 3x3 of the OpenGL style "view" transformation from the transpose of the camera orientation
	// This is the inversion process. Since these 3x3 matrices are orthonormal a transpose is sufficient to invert
//	m3dLoadIdentity44f(modelingTransform);	
//	for (int i = 0; i < 3; i++) {
//		for (int j = 0; j < 3; j++) {
//			MatrixElement(modelingTransform, i, j) = MatrixElement(openGLModelViewTransform, j, i);
//		}
//	}

	
	
	
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
	
	g_x_label.text	= [NSString stringWithFormat:@"%.2f", g[0]];
	g_y_label.text	= [NSString stringWithFormat:@"%.2f", g[1]];
	g_z_label.text	= [NSString stringWithFormat:@"%.2f", g[2]];
	
	
	

	
	
	
	
	// Update G
	m3dCopyVector3f(gPast, g);

	
}

@end
