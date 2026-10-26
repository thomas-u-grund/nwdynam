*! version 1.1.3  26oct2026 _nwdynam_libcheck - make the nwdynam Mata library visible
* Copyright (C) 2026 Thomas U. Grund
* This file is part of nwdynam (GPL-3, see the file LICENSE).
* Stata indexes the Mata libraries (.mlib) once, early in a session. A library installed by net install can be missing from that
* index (then Mata stops with r(3499) "<function>() not found"). nwdynam calls this program first: when a function of its library is
* not visible, it runs -mata mlib index- once.
capture program drop _nwdynam_libcheck
program define _nwdynam_libcheck
	version 14
	capture mata: mata which DynamFitUnit1()
	if _rc {
		capture mata: mata mlib index
	}
end
