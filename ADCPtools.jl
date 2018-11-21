module ADCPtools

export hello, getzaxmap

hello(a::AbstractString) = "Hello, $a"

# function getzaxmap(rbeam::Float32, theta::Float32)
function getzaxmap(rbeam, theta)
# USAGE
# -----
# mapped_zaxis = getzaxmap(rbeam, theta)
#
# Returns the projection of the along-beam coordinate 'rbeam'
# from the local vertical axis. The along-beam direction is tilted
# 'theta' degrees from the vertical.
  Zmap = rbeam.*cos(theta.*pi/180)
  dZu = Zmap[end] - Zmap[end-1]
  rmax = rbeam[end]
  nz = length(Zmap)
  while Zmap[end]<(rmax-dZu)
    Zmap = [Zmap Zmap[end]+dZu] # Complete vertical axis to map to.
  end
end

end
