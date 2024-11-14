module OffsetArraysExt
import DSP
import OffsetArrays

DSP.conv_axis_with_offset(::OffsetArrays.IdOffsetRange) = true

end
