//
//  HelloAccelerometerAppDelegate.h
//  HelloAccelerometer
//
//  Created by turner on 5/6/09.
//  Copyright Douglass Turner Consulting 2009. All rights reserved.
//

#import <UIKit/UIKit.h>

@class HelloAccelerometerController;

@interface HelloAccelerometerAppDelegate : NSObject <UIApplicationDelegate> {
	
	HelloAccelerometerController *controller;
    UIWindow *window;
}

@property (nonatomic, retain) IBOutlet UIWindow *window;
@property (nonatomic, retain) IBOutlet HelloAccelerometerController *controller;

@end

