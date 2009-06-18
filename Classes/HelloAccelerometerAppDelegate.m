//
//  HelloAccelerometerAppDelegate.m
//  HelloAccelerometer
//
//  Created by turner on 5/6/09.
//  Copyright Douglass Turner Consulting 2009. All rights reserved.
//

#import "HelloAccelerometerAppDelegate.h"
#import "HelloAccelerometerController.h"

@implementation HelloAccelerometerAppDelegate

@synthesize window;
@synthesize controller;

- (void)applicationDidFinishLaunching:(UIApplication *)application {    
	
	[controller initializeAccelerometer];
	
    [window makeKeyAndVisible];
}


- (void)dealloc {
    [window release];
    [super dealloc];
}


@end
