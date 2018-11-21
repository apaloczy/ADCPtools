function idxfl = nearfl(x, x0)
  % Get the lowest index of the two points in
  % vector x that bound the number x0.
  dx = abs(x - x0);
  flr = mink(dx, 2);
  il = find(dx==flr(1));
  ir = find(dx==flr(2));
  idxfl = min([il ir]);
end
