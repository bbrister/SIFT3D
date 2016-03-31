classdef Sift3DParser < inputParser
    %Sift3DParser helper class to parse SIFT3D options. This class inherits
    % from inputParser. Call its 'parseAndVerify' method to parse the input 
    % string. The results will be stored in the 'Results' field. Additional
    % options can be added just as in inputParser.
    
    properties (SetAccess = private)
        % Option names
        peakThreshStr = 'peakThresh';
        cornerThreshStr = 'cornerThresh';
        numKpLevelsStr = 'numKpLevels';
        sigmaNStr = 'sigmaN';
        sigma0Str = 'sigma0';
    end
    
    methods
        
        % Constructor
        function self = Sift3DParser
            
            % Call the parent constructor
            self = self@inputParser;
            
            % Add the SIFT3D options
            self.addParamValue(self.peakThreshStr, [])
            self.addParamValue(self.cornerThreshStr, [])
            self.addParamValue(self.numKpLevelsStr, [])
            self.addParamValue(self.sigmaNStr, [])
            self.addParamValue(self.sigma0Str, [])
        end
        
        % Verify and retrieve the SIFT3D options
        function optStruct = parseAndVerify(self, varargin)
            
            % Parse the input
            self.parse(varargin{:})
              
            % Retrieve the results
            peakThresh = self.Results.peakThresh;
            cornerThresh = self.Results.cornerThresh;
            numKpLevels = self.Results.numKpLevels;
            sigmaN = self.Results.sigmaN;
            sigma0 = self.Results.sigma0;
            
            % Verify the results
            if ~isempty(peakThresh)
                validateattributes(peakThresh, {'numeric'}, ...
                    {'real', 'positive', 'scalar', '<=', 1}, 'peakThresh')
            end
            if ~isempty(cornerThresh)
                validateattributes(cornerThresh, {'numeric'}, ...
                    {'real', 'nonnegative', 'scalar', '<=', 1}, ...
                    'cornerThresh')
            end
            if ~isempty(numKpLevels)
                validateattributes(numKpLevels, {'numeric'}, ...
                    {'real', 'integer', 'scalar', 'positive'}, ...
                    'numKpLevels')
            end
            if ~isempty(sigmaN)
                validateattributes(sigmaN, {'numeric'}, ...
                    {'real', 'positive', 'scalar'}, 'sigmaN')
            end
            if ~isempty(sigma0)
                validateattributes(sigma0, {'numeric'}, ...
                    {'real', 'positive', 'scalar'}, 'sigma0')
            end
            if sigmaN >= sigma0
                error('Cannot have sigmaN >= sigma0')
            end            
            
            % Collect the options in a struct
            optStruct = struct( ...
                self.peakThreshStr, peakThresh, ...
                self.cornerThreshStr, cornerThresh, ...
                self.numKpLevelsStr, numKpLevels, ...
                self.sigmaNStr, sigmaN, ...
                self.sigma0Str, sigma0);
        end
    end
end