* Common set-up of the example do-files of the manual. The environment variables
* NWCORE (directory of the nwcommands core), NWDYNAM (directory of nwdynam) and
* NWPLUS (an empty scratch directory used as PLUS) are set by run_examples.sh.
capture program drop nwinstall
local nwcore : env NWCORE
local nwdynam : env NWDYNAM
local nwplus : env NWPLUS
sysdir set PLUS "`nwplus'"
quietly adopath ++ "`nwcore'"
quietly adopath ++ "`nwcore'/lib"
quietly adopath ++ "`nwdynam'/lib"
quietly adopath ++ "`nwdynam'/lib/plugins/macos"
quietly adopath ++ "`nwdynam'"
set more off
set linesize 80

* mkmail #: declares the emails of the first # days of the manufacturing data as the event network mail,
* with every employee of the company in the risk set (the code is shown in full in ex01_first.do)
capture program drop mkmail
program mkmail
    args tmax
    quietly {
        nwwebuse manufacturing, nwclear
        keep if t < `tmax'
        levelsof sender, local(a) clean
        levelsof receiver, local(b) clean
        local seen : list a | b
        tempfile mail
        save `mail'
        nwwebuse manufacturing_nodes, nwclear
        levelsof label, local(everyone) clean
        local without : list everyone - seen
        local miss
        foreach e of local without {
            local miss `miss',`e'
        }
        local miss = substr("`miss'", 2, .)
        use `mail', clear
        nwset sender receiver, directed eventtime(t) name(mail)
        nwaddnodes mail, nodenames(`miss')
    }
end

