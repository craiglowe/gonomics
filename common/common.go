package common

import (
	"github.com/golang/glog"
)

func ExitIfError(err error) {
	if err != nil {
		glog.Exit(err)
	}
}

func Exit(message string) {
	glog.Exit(message)
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

func TripleMax(a int, b int, c int) int {
	if a >= b && a >= c {
		return a
	} else if b >= c {
		return b
	} else {
		return c
	}
}

func TripleMin(a int, b int, c int) int {
	if a <= b && a <= c {
		return a
	} else if b <= c {
		return b
	} else {
		return c
	}
}
