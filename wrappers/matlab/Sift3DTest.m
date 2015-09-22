classdef Sift3DTest < TestCase
    %Sift3DTest a test suite for SIFT3D.
    %
    % To run this test suite, you must have installed xUnit. As of August
    % 25th, 2015, xUnit is available at:
    %   http://www.mathworks.com/matlabcentral/fileexchange/22846-matlab-xunit-test-framework
    %
    % Run the tests with the following command:
    %   runtests
    %
    % This test suite can only be run from the build tree.
    %
    % Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
    
    properties (SetAccess = private)
        cd
        buildDir
        binDir
        examplesDir
        im1Name
        im2Name
        dataName
        kpCmd
        tolText
    end
    
    methods
        
        % Constructor
        function self = Sift3DTest(name)
            
            % Call the parent constructor
            self = self@TestCase(name);
            
            % Current directory
            self.cd = fileparts(mfilename('fullpath'));
            
            % Build directory
            self.buildDir = fullfile(self.cd, '..', '..', '..');
            
            % Binary directory
            self.binDir = fullfile(self.buildDir, 'bin');
            
            % Examples directory
            self.examplesDir = fullfile(self.buildDir, 'examples');
            
            % Image file names
            self.im1Name = fullfile(self.examplesDir, '1.nii.gz');
            self.im2Name = fullfile(self.examplesDir, '2.nii.gz');
            
            % Data file name
            self.dataName = fullfile(self.examplesDir, 'data.mat');
            
            % Keypoints command name
            self.kpCmd = fullfile(self.binDir, 'kpSift3D');
            
            % Error tolerance for text output
            self.tolText = 0.01;
            
        end
        
        %         % Test keypoint detection against the CLI version
        %         function detectCliTest(self)
        %
        %             % Output file name
        %             kpCliName = 'kpCli.csv';
        %
        %             % Detect keypoints using the command line interface
        %             status = runCmd([self.kpCmd ' --keys ' kpCliName ' ' ...
        %                 self.im1Name]);
        %             assertEqual(status, 0);
        %
        %             % Load the CLI keypoints
        %             kpCli = csvread(kpCliName);
        %
        %             % Load the image data
        %             load(self.dataName);
        %
        %             % Detect keypoints using matlab
        %             keys = detectSift3D(im1);
        %
        %             % Check the dimensions
        %             assertEqual(size(kpCli, 1), length(keys));
        %             assertEqual(size(kpCli, 2), numel(keys(1).coords) + ...
        %                 numel(keys(1).scale) + numel(keys(1).ori));
        %
        %             % Compare the two
        %             for i = 1 : length(keys)
        %
        %                 mKey = keys(i);
        %                 cliKey = kpCli(i, :);
        %
        %                 % Check the coordinates
        %                 assertElementsAlmostEqual(mKey.coords, cliKey(1:3), ...
        %                     'absolute', self.tolText);
        %
        %                 % Check the scale
        %                 assertElementsAlmostEqual(mKey.scale, cliKey(4), ...
        %                     'absolute', self.tolText);
        %
        %                 % Check the orientation
        %                 assertElementsAlmostEqual(mKey.ori, ...
        %                     reshape(cliKey(5:end), size(mKey.ori))', ...
        %                     'absolute', self.tolText);
        %             end
        %
        %             % Clean up
        %             delete(kpCliName);
        %         end
        %
        %         % Test descriptor extraction against the CLI version
        %         function extractCliTest(self)
        %
        %             % Output file name
        %             descCliName = 'descCli.csv';
        %
        %             % Extract descriptors using the command line interface
        %             status = runCmd([self.kpCmd ' --desc ' descCliName ' ' ...
        %                 self.im1Name]);
        %             assertEqual(status, 0);
        %
        %             % Read the results
        %             descCli = csvread(descCliName);
        %
        %             % Load the image data
        %             load(self.dataName);
        %
        %             % Extract descriptors using matlab
        %             keys = detectSift3D(im1);
        %             [desc, coords] = extractSift3D(keys);
        %
        %             % Check the dimensions
        %             assertEqual(size(desc, 1), size(coords, 1));
        %             assertEqual(size(descCli, 1), size(desc, 1));
        %             assertEqual(size(descCli, 2), size(desc, 2) + size(coords, 2));
        %
        %             % Compare the two
        %             for i = 1 : length(keys)
        %
        %                 cliDescrip = descCli(i, :);
        %
        %                 % Check the coordinates
        %                 assertElementsAlmostEqual(cliDescrip(1 : 3), ...
        %                     coords(i, :), 'absolute', self.tolText);
        %
        %                 % Check the descriptor
        %                 assertElementsAlmostEqual(cliDescrip(4 : end), ...
        %                     desc(i, :), 'absolute', self.tolText);
        %             end
        %
        %             % Clean up
        %             delete(descCliName);
        %         end
        %
        %         % Test that "raw" image descriptors are close to those extracted
        %         % from a Gaussian scale-space pyramid
        %         function rawTest(self)
        %
        %             % Load the image data
        %             load(self.dataName);
        %
        %             % Detect keypoints
        %             keys = detectSift3D(im1);
        %
        %             % Extract descriptors using the pyramid
        %             [descPyr, coordsPyr] = extractSift3D(keys);
        %
        %             % Extract raw descriptors
        %             [descRaw, coordsRaw] = extractSift3D(keys, im1);
        %
        %             % Check the results
        %             assertElementsAlmostEqual(coordsPyr, coordsRaw);
        %             assertElementsAlmostEqual(descPyr, descRaw, 'absolute', 0.2);
        %
        %         end
        
        % Test reading and writing a NIFTI image
        function niftiIOTest(self)
            
            % The temporary file name
            imName = 'temp.nii.gz';
            
            % Make random image data
            imWritten = rand(10, 15, 20);
            
            % Write the image as a NIFTI file
            imWrite3D(imName, imWritten);
            
            % Read the image back
            imRead = imRead3D(imName);
            
            % Clean up
            delete(imName);
            
            % Ensure the results are identical
            assertElementsAlmostEqual(imWritten, imRead, 'relative', 1E-3);
        end
        
        % Test reading and writing a DICOM image
        function dicomIOTest(self)
            
            % The temporary file name
            imName = 'temp.dcm';
            
            % Make random image data
            imWritten = rand(10, 15, 20);
            
            % Write the image as a DICOM file
            imWrite3D(imName, imWritten);
            
            % Read the image back
            imRead = imRead3D(imName);
            
            % Clean up
            delete(imName);
            
            % Ensure the results are identical
            assertElementsAlmostEqual(imWritten, imRead, 'relative', 1E-3);
        end
        
        % Test reading and writing a directory of DICOM images
        function dirIOTest(self)
            
            % The temporary file name
            dirName = 'temp';
            
            % Make random image data
            imWritten = rand(10, 15, 20);
            
            % Write the image as a DICOM file
            imWrite3D(dirName, imWritten);
            
            % Read the image back
            imRead = imRead3D(dirName);
            
            % Clean up
            delete(dirName);
            
            % Ensure the results are identical
            assertElementsAlmostEqual(imWritten, imRead, 'relative', 1E-3);
        end
        
        % Test reading and writing a 2D image
        function io2DTest(self)
             % The temporary file name
            imName = 'temp.nii.gz';
            
            % Make random image data
            imWritten = rand(10, 15);
            
            % Write the image as a NIFTI file
            imWrite3D(imName, imWritten);
            
            % Read the image back
            imRead = imRead3D(imName);
            
            % Clean up
            delete(imName);
            
            % Ensure the results are identical
            assertElementsAlmostEqual(imWritten, imRead, 'relative', 1E-3);
        end
        
        % Test reading an invalid filetype
        function readInvalidTypeTest(self)
            
            % Temporary file name
            tempFileName = 'temp.mat';
            
            % Make a temporary file
            threwErr = false;
            save(tempFileName);
            
            % Try to read it
            try
                im = imRead3D(tempFileName);
            catch ME
                threwErr = true;
            end
            
            % Clean up
            delete(tempFileName);
            
            assertTrue(threwErr);
        end
        
        % Test writing an invalid filetype
        function writeInvalidTypeTest(self)
            im = rand(10, 10, 10, 10);
            threwErr = false;
            try
                imWrite3D(im, 'bogus.txt')
            catch ME
                threwErr = true;
            end
            assertTrue(threwErr);
        end
        
        % Test reading from a nonexistent file
        function imReadNonexistentTest(self)
            
            threwErr = false;
            try
                im = imRead3D('nonexistent.nii.gz');
            catch ME
                threwErr = true;
            end
            assertTrue(threwErr);
        end
    end
end

function status = runCmd(cmd)
%runCmd helper function to run a command in the default system environment,
% without Matlab's changes to the environment variables

% Strip the LD_LIBRARY_PATH environment variable of Matlab
% directories
ldPath = 'LD_LIBRARY_PATH';
oldLdPath = getenv(ldPath);
newLdPath = regexprep(oldLdPath, '[^:]*MATLAB[^:]*:*', '');
setenv(ldPath, newLdPath);

% Run the command
status = system(cmd);

% Return the environment to its previous state
setenv(ldPath, oldLdPath);

end
