# ================================================
# VMD Atom Renamer with One Persistent Dialog
# ================================================

# 1) Prepare your molecule
set molid      0
set num_atoms [molinfo $molid get numatoms]
set atom_index 0

# 2) Create the two representations
# mol delrep 0 $molid
# mol representation Lines
# mol selection "all"
# mol addrep $molid

mol representation VDW 0.1
mol color Name
mol selection "index 0"
mol addrep $molid

# 3) Build a persistent dialog window
toplevel .dlg
wm title .dlg "Atom Renamer"

# Info label
label .dlg.info -justify left -anchor w -width 30 -text ""
pack .dlg.info -padx 10 -pady 5

# Entry for new name
entry .dlg.entry -width 10
pack .dlg.entry -padx 10
bind .dlg.entry <Return> { perform_rename }

# Buttons
frame .dlg.btns
button .dlg.btns.next  -text "Next (Enter)" -command perform_rename
button .dlg.btns.exit  -text "Exit"         -command { destroy .dlg; mol modselect 1 $molid "none"; puts "Exited at atom $::atom_index" }
pack .dlg.btns.next .dlg.btns.exit -side left -padx 5 -pady 5

# 4) Procedure to update dialog for current atom
proc update_dialog {} {
    global molid num_atoms atom_index

    if { $atom_index >= $num_atoms } {
        mol modselect 1 $molid "none"
        puts "\n✅ All atoms done."

        # ------------------------------------------------------------------
        # NEW: back‑up the original file, then overwrite it with new names
        # ------------------------------------------------------------------
        set origfile [molinfo $molid get filename]
        if { $origfile ne "" && [file exists $origfile] } {
            # choose a backup name that never clobbers an existing file
            set backup $origfile.bak
            if { [file exists $backup] } {
                set ts [clock format [clock seconds] -format "%Y%m%d%H%M%S"]
                set backup "${origfile}.${ts}.bak"
            }
            file copy -force $origfile $backup
            puts "📁  Original backed up to '$backup'"

            # write the renamed structure back over the original
            set selall [atomselect $molid "all"]
            $selall writepdb $origfile
            $selall delete
            puts "💾  Renamed structure saved as '$origfile'"
        } else {
            puts "⚠️  Molecule has no readable filename; nothing written."
        }
        # ------------------------------------------------------------------

        destroy .dlg
        return
    }

    # Highlight current atom
    mol modselect 1 $molid "index $atom_index"
    update

    # Get atom info
    set sel      [atomselect $molid "index $atom_index"]
    set oldname  [lindex [$sel get name]    0]
    set resid    [lindex [$sel get resid]   0]
    set resname  [lindex [$sel get resname] 0]
    set element  [lindex [$sel get element] 0]
    $sel delete

    # Update label and entry
    .dlg.info configure -text "Atom #$atom_index of $num_atoms\nResid $resid ($resname)\nElement $element\nCurrent name: $oldname"
    .dlg.entry delete 0 end
    .dlg.entry insert 0 $oldname
    focus .dlg.entry
}

# 5) Procedure called when Enter or Next is pressed
proc perform_rename {} {
    global molid atom_index

    # read new name
    set newname [.dlg.entry get]
    # get old name again for comparison
    set sel   [atomselect $molid "index $atom_index"]
    set oldnm [lindex [$sel get name] 0]
    if { $newname ne $oldnm } {
        $sel set name $newname
        puts "Atom $atom_index renamed to '$newname'"
    } else {
        puts "Atom $atom_index left as '$oldnm'"
    }
    $sel delete

    incr atom_index
    update_dialog
}

# 6) Start
update_dialog
