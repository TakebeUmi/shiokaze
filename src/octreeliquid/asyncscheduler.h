/*
**	asyncscheduler.h
**
**	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
**	Created by Ryoichi Ando <rand@nii.ac.jp> on November 30, 2019. All rights reserved.
**	Permission is hereby granted, free of charge, to any person obtaining a copy of
**	this software and associated documentation files (the "Software"), to deal in
**	the Software without restriction, including without limitation the rights to use,
**	copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the
**	Software, and to permit persons to whom the Software is furnished to do so,
**	subject to the following conditions:
**
**	The above copyright notice and this permission notice shall be included in all copies
**	or substantial portions of the Software.
**
**	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
**	INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
**	PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
**	HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
**	CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
**	OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
//
#ifndef SHKZ_SCHEDULER_H
#define SHKZ_SCHEDULER_H
//
#include <shiokaze/core/configurable.h>
//
SHKZ_BEGIN_NAMESPACE
//
class asyncscheduler : public recursive_configurable, public credit {
private:
	//
	LONG_NAME("Asynchronous Scheduler")
	ARGUMENT_NAME("AsyncScheduler")
	//
	virtual void configure( configuration &config ) override {
	}
	//
	asyncscheduler ( recursive_configurable *parent ) {
		if( parent ) parent->add_child(this);
		else setup_now();
	}
};
//
SHKZ_END_NAMESPACE
//
#endif
//