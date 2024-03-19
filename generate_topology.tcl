package require autopsf

# Define the main directory
set mainDir "./CASP13_BB"

# Create the output directory for psfgen files
set psfgenDir "${mainDir}_PSF"
file mkdir $psfgenDir
# Get a list of subfolders in the main directory
set subfolders [glob -type d $mainDir/*]

set topology top_all36_prot.rtf
# Loop through each subfolder
foreach subfolder $subfolders {
    # Extract the folder name from the path
    set foldername [file tail $subfolder]

    # Specify the input PDB file for each subfolder
    set pdbFile "${subfolder}/1bb1.pdb"

    # Specify the output PSF file for each subfolder
    set prefix "${psfgenDir}/${foldername}/"

    # Create the output directory for each subfolder
    file mkdir [file join $psfgenDir $foldername]
    set id [mol new $pdbFile]
    mol rename $id tmp
    autopsf -mol $id -prefix tmp -top $topology
    file rename tmp_formatted_autopsf.pdb ${prefix}/1psf.pdb
    file delete tmp_formatted_autopsf.psf tmp_formatted_autopsf.log tmp_formatted.pdb
    mol delete $id

}

puts "PSF generation completed."
exit