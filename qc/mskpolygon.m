function varargout = mskpolygon(x, y, A)
  % USAGE
  % -----
  % Aq = mskpolygon(x, y, A)
  %
  %       --OR--
  %
  % [Aq, msks] = mskpolygon(x, y, A)
  %
  % Given an array `A(x, y)`, select polygon vertices using points selected
  % using ginput, and return a copy of `A` (`Aq`) with all points inside
  % any one of the polygon(s) NaN'd out. `msks` is a 3D stack of logical
  % arrays, where each page indicates one of the polygons selected.
  %
  % A *left click* selects the vertices of an arbitrary polygon. Placing the cursor
  % on the chosen spot and pressing one of the h' or 'v' keys on the keyboard
  % indicates the following actions:
  %
  % Keyboard/mouse char/button  |   action
  %-----------------------------------------------------------------------------
  %     'left click'            | define vertices of arbitrary polygon
  %         'h'                 | define lateral edges of _vertical_ rectangle
  %         'v'                 | define vertical edges of _horizontal_ rectangle
  %
  %
  % *Right-click* to indicate the end of the points defining each polygon,
  % and start drawing the next one. Press *enter* to finish.
  %
  % Author: André Palóczy (paloczy@gmail.com)
  % Nov/23/2018
  if nargout>2
    error('Invalid number of output arguments')
  end

  if length(x)==numel(x) && length(y)==numel(y) % 1D case.
    [x, y] = meshgrid(x, y);
  elseif ndims(x)==2 && ndims(y)==2
    ;
  else
    error('Inconsistent shape in `x` and/or `y` arrays')
  end

  [xx, yy, wat] = ginput;
  wat = [wat; 3];

  delim = find(wat==3);
  npoly = numel(delim);

  Aq = A;
  msks = [];
  for n=1:npoly
    d = delim(n);
    w = wat(1:d-1);
    fa = find(w==1);   % 'left-click'
    fh = find(w==104); % 'h'
    fv = find(w==118); % 'v'
    xa = xx(fa); ya = yy(fa);
    xh = sort(xx(fh)); yh = sort(yy(fh));
    xv = sort(xx(fv)); yv = sort(yy(fv));

    if numel(fa)>=3
      xx = [xx; xx(1)]; yy = [yy; yy(1)]; w = [w; w(1)]; % Close polygon.
      in = inpolygon(x, y, xa, ya);
      Aq(in) = NaN;
      if nargout==2 % keep masks.
        if isempty(msks)
          msks = in;
        else
          msks = cat(3, msks, in);
        end
      end
    elseif isempty(fa)
      ;
    else % There should be at least 3 points of kind 'a':
      error('There must be at least 3 points for arbitrary rectangles.')
    end

    if numel(fh)==2
      in = x>=xh(1) & x<=xh(2);
      Aq(in) = NaN;
      if nargout==2 % keep maskses
        if isempty(msks)
          msks = in;
        else
          msks = cat(3, msks, in);
        end
      end
    elseif isempty(fh)
      ;
    else % There should be only two points of kind 'h' and 'v':
      error('There must be exactly 2 points for the horizontal and vertical rectangles.')
    end

    if numel(fv)==2
      in = y>=yv(1) & y<=yv(2);
      Aq(in) = NaN;
      if nargout==2 % keep maskses
        if isempty(msks)
          msks = in;
        else
          msks = cat(3, msks, in);
        end
      end
    elseif isempty(fv)
      ;
    else % There should be only two points of kind 'h' and 'v':
      error('There must be exactly 2 points for the horizontal and vertical rectangles.')
    end

    wat = wat(d+1:end);
    delim = delim - d;
  end % for

  varargout{1} = Aq;
  if nargout==2
    varargout{2} = msks;
  end

end
