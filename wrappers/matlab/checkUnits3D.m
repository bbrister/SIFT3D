function units = checkUnits3D(units)
%checkUnits3D helper routine to check the validity of the 'units' argument
% to various SIFT3D functions. Returns a properly formatted version of 
% units.

if ~isvector(units) || length(units) ~= 3 || ~isnumeric(units) || ...
    ~isreal(units)
    error('units must be a [3x1] real numeric vector')
end

if any(units <= 0)
    error('units must be positive')
end

if isequal(size(units), [1 3])
    units = units';
end

units = double(units);

end
