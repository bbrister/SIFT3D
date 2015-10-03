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
        fullTest
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
            
            % Keypoints command name
            self.kpCmd = fullfile(self.binDir, 'kpSift3D');
            
            % Error tolerance for text output
            self.tolText = 0.01;
            
            % Run the tests on real data (slow)
            fullTest = true;
            
        end
        
        % Test keypoint detection against the CLI version
        function detectCliTest(self)
            
            if (self.fullTest)
                
                % Output file name
                kpCliName = 'kpCli.csv';
                
                % Detect keypoints using the command line interface
                status = runCmd([self.kpCmd ' --keys ' kpCliName ' ' ...
                    self.im1Name]);
                assertEqual(status, 0);
                
                % Load the CLI keypoints
                kpCli = csvread(kpCliName);
                
                % Load the image data
                im1 = imRead3D(self.im1Name);
                
                % Detect keypoints using matlab
                keys = detectSift3D(im1);
                
                % Check the dimensions
                assertEqual(size(kpCli, 1), length(keys));
                assertEqual(size(kpCli, 2), numel(keys(1).coords) + ...
                    numel(keys(1).scale) + numel(keys(1).ori));
                
                % Compare the two
                for i = 1 : length(keys)
                    
                    mKey = keys(i);
                    cliKey = kpCli(i, :);
                    
                    % Check the coordinates
                    assertElementsAlmostEqual(mKey.coords, cliKey(1:3), ...
                        'absolute', self.tolText);
                    
                    % Check the scale
                    assertElementsAlmostEqual(mKey.scale, cliKey(4), ...
                        'absolute', self.tolText);
                    
                    % Check the orientation
                    assertElementsAlmostEqual(mKey.ori, ...
                        reshape(cliKey(5:end), size(mKey.ori))', ...
                        'absolute', self.tolText);
                end
                
                % Clean up
                delete(kpCliName);
                
            end
        end
        
        % Test descriptor extraction against the CLI version
        function extractCliTest(self)
            
            if (self.fullTest)
                
                % Output file name
                descCliName = 'descCli.csv';
                
                % Extract descriptors using the command line interface
                status = runCmd([self.kpCmd ' --desc ' descCliName ' ' ...
                    self.im1Name]);
                assertEqual(status, 0);
                
                % Read the results
                descCli = csvread(descCliName);
                
                % Load the image data
                im1 = imRead3D(self.im1Name);
                
                % Extract descriptors using matlab
                keys = detectSift3D(im1);
                [desc, coords] = extractSift3D(keys);
                
                % Check the dimensions
                assertEqual(size(desc, 1), size(coords, 1));
                assertEqual(size(descCli, 1), size(desc, 1));
                assertEqual(size(descCli, 2), size(desc, 2) + size(coords, 2));
                
                % Compare the two
                for i = 1 : length(keys)
                    
                    cliDescrip = descCli(i, :);
                    
                    % Check the coordinates
                    assertElementsAlmostEqual(cliDescrip(1 : 3), ...
                        coords(i, :), 'absolute', self.tolText);
                    
                    % Check the descriptor
                    assertElementsAlmostEqual(cliDescrip(4 : end), ...
                        desc(i, :), 'absolute', self.tolText);
                end
                
                % Clean up
                delete(descCliName);
                
            end
            
        end
        
        % Test that "raw" image descriptors are close to those extracted
        % from a Gaussian scale-space pyramid
        function rawTest(self)
            
            if (self.fullTest)
                
                % Load the image data
                im1 = imRead3D(self.im1Name);
                
                % Detect keypoints
                keys = detectSift3D(im1);
                
                % Extract descriptors using the pyramid
                [descPyr, coordsPyr] = extractSift3D(keys);
                
                % Extract raw descriptors
                [descRaw, coordsRaw] = extractSift3D(keys, im1);
                
                % Check the results
                assertElementsAlmostEqual(coordsPyr, coordsRaw);
                assertElementsAlmostEqual(descPyr, descRaw, 'absolute', 0.2);
                
            end
            
        end
        
        % Test extracting descriptors with no input image, before keypoints
        % were detected
        function extractBeforeDetectTest(self)
            
            keys = keypoint3D([1 1 1]);
            
            threwErr = false;
            try
                desc = extractSift3D(keys);
            catch ME
                threwErr = true;
            end
            assertTrue(threwErr);
        end
        
        % Test reading and writing a NIFTI image
        function niftiIOTest(self)
            
            % The temporary file name
            imName = 'temp.nii.gz';
            
            % Make random image data, scaled to the range [0, 1]
            imWritten = rand(10, 15, 20);
            imWritten = imWritten / max(imWritten(:));
            
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
            
            % Remove any past instances of this file
            if exist(imName, 'file')
                delete(imName)
            end
            
            % Make random image data, scaled to the range [0, 1]
            imWritten = rand(10, 15, 20);
            imWritten = imWritten / max(imWritten(:));
            
            % Write the image as a DICOM file
            imWrite3D(imName, imWritten);
            
            % Read the image back
            imRead = imRead3D(imName);
            
            % Clean up
            delete(imName);
            
            % Ensure the results are identical
            assertElementsAlmostEqual(imWritten, imRead, 'absolute', 1E-2);
        end
        
        % Test reading and writing a directory of DICOM images
        function dirIOTest(self)
            
            % The temporary file name
            dirName = 'temp';
            
            % Make random image data, scaled to the range [0, 1]
            imWritten = rand(10, 15, 20);
            imWritten = imWritten / max(imWritten(:));
            
            % Write the image as a DICOM file
            imWrite3D(dirName, imWritten);
            
            % Read the image back
            imRead = imRead3D(dirName);
            
            % Clean up
            rmdir(dirName, 's');
            
            % Ensure the results are identical
            assertElementsAlmostEqual(imWritten, imRead, 'absolute', 1E-2);
        end
        
        % Test reading and writing a 2D image
        function io2DTest(self)
            % The temporary file name
            imName = 'temp.nii.gz';
            
            % Make random image data, scaled to the range [0, 1]
            imWritten = rand(10, 15, 20);
            imWritten = imWritten / max(imWritten(:));
            
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
        
        % Test making valid keypoints
        function keypointValidTest(self)
            % Make the keypoints
            coords = [1 1 1; 2 2 2];
            scale = [1 2];
            ori = repmat(rotMat(1), [1 1 length(scale)]);
            keys = keypoint3D(coords, scale, ori);
            
            % Test the results
            for i = 1 : length(keys)
                key = keys(i);
                assertEqual(key.coords, coords(i, :));
                assertEqual(key.scale, scale(i));
                assertEqual(key.ori, ori(:, :, i));
            end
        end
        
        % Test making keypoints with invalid coordinate dimensions
        function keypointInvalidCoordTest(self)
            % Make the keypoints
            coords = [1 1 1 1];
            
            % Test the results
            threwErr = false;
            try
                keypoint3D(coords);
            catch ME
                threwErr = true;
            end
            assertTrue(threwErr);
        end
        
        % Test making keypoints with invalid scale dimensions
        function keypointInvalidScaleTest(self)
            % Make the keypoints
            coords = [1 1 1; 2 2 2];
            scale = [1 1; 2 2];
            
            % Test the results
            threwErr = false;
            try
                keypoint3D(coords, scale);
            catch ME
                threwErr = true;
            end
            assertTrue(threwErr);
        end
        
        % Test making keypoints with invalid rotation matrix dimensions
        function keypointInvalidRotDimsTest(self)
            % Make the keypoints
            coords = [1 1 1; 2 2 2];
            ori = repmat(ones(3, 4), [1 1 size(coords, 1)]);
            
            % Test the results
            threwErr = false;
            try
                keypoint3D(coords, [], ori);
            catch ME
                threwErr = true;
            end
            assertTrue(threwErr);
        end
        
        % Test making a keypoint with a reflection matrix
        function keypointReflectTest(self)
            % Make the keypoints
            coords = [1 1 1];
            ori = eye(3);
            ori(1, 1) = -1;
            
            % Test the results
            threwErr = false;
            try
                keypoint3D(coords, [], ori);
            catch ME
                threwErr = true;
            end
            assertTrue(threwErr);
        end
        
        % Test making a keypoint with a non-orthogonal matrix
        function keypointOrthTest(self)
            % Make the keypoints
            coords = [1 1 1];
            ori = rand(3);
            
            % Ensure that the determinant is equal to one, so it is
            % orthogonality and not reflection which causes an error
            detOri = det(ori);
            factor = sign(detOri) * abs(detOri) ^ (-1 / length(ori));
            ori = ori * factor;
            assert(abs(det(ori) - 1) < eps * 1E2);
            
            % Test the results
            threwErr = false;
            try
                keypoint3D(coords, [], ori);
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
ldPathVar = 'LD_LIBRARY_PATH';
oldLdPath = getenv(ldPathVar);
newLdPath = regexprep(oldLdPath, '[^:]*MATLAB[^:]*:*', '');
setenv(ldPathVar, newLdPath);

% Run the command
status = system(cmd);

% Return the environment to its previous state
setenv(ldPathVar, oldLdPath);

end

function R = rotMat(theta)
%rotMat Helper function to make a 3D rotation matrix, rotating by angle
% theta in the XY plane
R = eye(3);
R(1, 1) = cos(theta);
R(1, 2) = -sin(theta);
R(2, 1) = sin(theta);
R(2, 2) = cos(theta);
end
