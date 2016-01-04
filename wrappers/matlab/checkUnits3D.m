function units = checkUnits3D(units, name)
%checkUnits3D helper function to check the validity of the 'units' argument
% to various SIFT3D functions. Returns a properly formatted version of 
% units, with missing values set to 1.
%
% The optional 'name' argument is used in error messages.

% Default arguments
if nargin < 1 || isempty(units)
    units = ones(3, 1);
end
if nargin < 2 || isempty(name)
    name = 'units';
end

% Number of image dimensions
ndim = 3;

if ~isvector(units) || length(units) > ndim || ~isnumeric(units) || ...
    ~isreal(units)
    error([name ' must be a [3x1] real numeric vector'])
end

if any(units <= 0)
    error([name ' must be positive'])
end

% Transpose column vectors
if size(units, 2) > 1
    units = units';
end

% Convert to double
units = double(units);

% Default for missing values
if length(units) < ndim
    units(length(units) + 1 : ndim) = 1;
end

end