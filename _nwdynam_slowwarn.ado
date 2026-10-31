*! version 1.2.1  31oct2026 _nwdynam_slowwarn: warning line when the Mata fallback replaces the compiled plugin
* Copyright (C) 2026 Thomas U. Grund
* This file is part of nwdynam, a Stata package derived from the R package 'goldfish' (James Hollway, Christoph Stadtfeld, Marion Hoffman, Alvaro Uzaheta and contributors; GPL (>= 3)).
* nwdynam is free software: you can redistribute it and/or modify it under the terms of the GNU General
* Public License, version 3, as published by the Free Software Foundation. It is distributed in the hope
* that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
* or FITNESS FOR A PARTICULAR PURPOSE. See the file LICENSE for the full text.
* _nwdynam_slowwarn: one warning line (once per nwdynam call) before a fit whose engine could use the compiled plugin
* but the plugin file is missing for the platform, so the slower Mata engine is used (same wording style as nwergm).
capture program drop _nwdynam_slowwarn
program _nwdynam_slowwarn
	local os = cond(c(os) == "Windows", "windows", cond(c(os) == "MacOSX", "macos", "unix"))
	di as txt "(Warning: the compiled plugin is not used for this fit (the plugin file dynam_sim_`os'.plugin was not found; reinstall the package with {bf:net install nwdynam, replace}); the Mata fallback gives the same results but is slower. See {help nwdynam##native:Performance}.)"
end
