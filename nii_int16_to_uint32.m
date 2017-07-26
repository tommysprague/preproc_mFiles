% nii_int16_to_int32
%
% utility function for converting a 16-bit wraparound nii into a uint32 (as
% intended by Siemens)
%
% TCS 6/5/2017


function nii_int16_to_uint32(in_fn,out_fn)

mynii = niftiRead(in_fn);

mynii.data = mod(uint32(mynii.data),uint32(intmax('uint16')));

niftiWrite(mynii,out_fn);


return