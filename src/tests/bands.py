from abipy.abilab import abiopen
import abipy.data as abidata

# Open the file with energies computed on a k-path in the BZ
# and extract the band structure object.
with abiopen(abidata.ref_file("/home/wjing/eos_workflow/src/tests/job_2024-10-23-11-25-16-406622-41569/outdata/out_GSR.nc")) as nscf_file:
    nscf_ebands = nscf_file.ebands

# Open the file with energies computed with a homogeneous sampling of the BZ
# and extract the band structure object.
with abiopen(abidata.ref_file("/home/wjing/eos_workflow/src/tests/job_2024-10-23-11-22-18-082271-95854/outdata/out_GSR.nc")) as gs_file:
    gs_ebands = gs_file.ebands

edos = gs_ebands.get_edos()
nscf_ebands.plot_with_edos(edos, e0=None, with_gaps=True, title="Li XO Electron bands + DOS")
