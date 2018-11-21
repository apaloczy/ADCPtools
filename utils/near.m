function idx = near(x, x0)
  % Get the index of the point in
  % vector x nearest point closest
  % to the number x0.
  dx = abs(x - x0);
  fn = min(dx);
  idx = find(dx==fn);
end
