package common

import (
	"github.com/golang/glog"
)

func FatalIfError(err error) {
	if err != nil {
		glog.Exit(err)
	}
}

// zero is silent
// one is basic
// two is detailed
func Info(level glog.Level, s string) {
	glog.V(level).Info(s)
}

func Max(a int, b int) int {
	if a >= b {
		return a
	} else {
		return b
	}
}

func Min(a int, b int) int {
	if a <= b {
		return a
	} else {
		return b
	}
}
