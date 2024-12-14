module OffsetArraysExt
import DSP
import OffsetArrays

DSP.Convolutions.conv_axis_with_offset(::OffsetArrays.IdOffsetRange) = true

end
