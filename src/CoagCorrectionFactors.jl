module CoagCorrectionFactors

export intramodal_correction

"""
intramodal_correction(stdev)
Given the geometric standard deviation of an aerosol mode,
returns the correction factor for the intramodal coagulation integral
Correction factors for 0-th and 6-th moment are the same in Whitby 91 - see p. H.8
These correction factors are obtained from the CAM5 code for higher precision.
"""
function intramodal_correction(stdev)
      index = max(1, min(10, round(Integer, (stdev - 0.75)/0.25)))
      correction_factors = [
      0.707106785165097, 0.726148960080488, 0.766430744110958,
      0.814106389441342, 0.861679526483207, 0.903600509090092,   
      0.936578814219156, 0.960098926735545, 0.975646823342881,   
      0.985397173215326
      ]
      return correction_factors[index]
end

end
