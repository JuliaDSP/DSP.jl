module OffsetArraysExt
import DSP.Convolutions
import OffsetArrays

Convolutions.conv_axis_with_offset(::OffsetArrays.IdOffsetRange) = true

end
