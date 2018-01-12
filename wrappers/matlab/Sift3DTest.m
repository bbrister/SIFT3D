classdef Sift3DTest < TestCase
%Sift3DTest a test suite for SIFT3D.
%
% To run this test suite, you must install xUnit. As of January
% 12th, 2017, xUnit is available at:
% http://www.mathworks.com/matlabcentral/fileexchange/47302-xunit4
%
% Run the tests with the following command:
%   runxunit
%
% This test suite can only be run from the build tree.
%
% Copyright (c) 2015-2017 Blaine Rister et al., see LICENSE for details.

    properties (SetAccess = private)
        cd
        buildDir
        binDir
        examplesDir
        im1Name
        im2Name
        dataName
        kpCmd
        regCmd
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
            if ispc
                self.examplesDir = self.binDir;
            else
                self.examplesDir = fullfile(self.buildDir, 'examples');
            end
            
            % Image file names
            self.im1Name = fullfile(self.examplesDir, '1.nii.gz');
            self.im2Name = fullfile(self.examplesDir, '2.nii.gz');
            
            % Keypoints command name
            self.kpCmd = fullfile(self.binDir, 'kpSift3D');
            
            % Registration command name
            self.regCmd = fullfile(self.binDir, 'regSift3D');
            
            % Error tolerance for text output
            self.tolText = 0.01;
            
            % Run the tests on real data (slow)
            self.fullTest = true;
            
        end
        
        % Test keypoint detection against the CLI version
        function detectCliTest(self)
            
            if ~self.fullTest || ispc
                return
            end
                
            % Output file name
            kpCliName = 'kpCli.csv';
            
            % Detect keypoints using the command line interface
            status = runCmd([self.kpCmd ' --keys ' kpCliName ' ' ...
                self.im1Name]);
            assertEqual(status, 0);
            
            % Load the CLI keypoints
            kpCli = csvread(kpCliName);
            
            % Load the image data
            [im1, units1] = imRead3D(self.im1Name);
            
            % Detect keypoints using matlab
            keys = detectSift3D(im1, 'units', units1);
            
            % Check the dimensions
            assertEqual(size(kpCli, 1), length(keys));
            assertEqual(size(kpCli, 2), numel(keys(1).coords) + ...
                numel(keys(1).octave) + numel(keys(1).scale) + ...
                numel(keys(1).ori));
            
            % Compare the two
            for i = 1 : length(keys)
                
                mKey = keys(i);
                cliKey = kpCli(i, :);
                
                % Check the coordinates
                assertElementsAlmostEqual(mKey.coords, cliKey(1:3), ...
                    'absolute', self.tolText);
                
                % Check the octave
                assertEqual(mKey.octave, cliKey(4));
                
                % Check the scale
                assertElementsAlmostEqual(mKey.scale, cliKey(5), ...
                    'absolute', self.tolText);
                
                % Check the orientation
                assertElementsAlmostEqual(mKey.ori, ...
                    reshape(cliKey(6:end), size(mKey.ori))', ...
                    'absolute', self.tolText);
            end
            
            % Clean up
            delete(kpCliName);
        end
        
        % Test descriptor extraction against the CLI version
        function extractCliTest(self)
            
            if ~self.fullTest || ispc
                return
            end
                
                % Output file name
                descCliName = 'descCli.csv';
                
                % Extract descriptors using the command line interface
                status = runCmd([self.kpCmd ' --desc ' descCliName ' ' ...
                    self.im1Name]);
                assertEqual(status, 0);
                
                % Read the results
                descCli = csvread(descCliName);
                
                % Load the image data
                [im1, units1] = imRead3D(self.im1Name);
                
                % Extract descriptors using matlab
                keys = detectSift3D(im1, 'units', units1);
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
        
        % Test that "raw" image descriptors are close to those extracted
        % from a Gaussian scale-space pyramid
        function rawDescriptorTest(self)
            
            if ~self.fullTest
                return
            end
                
            % Load the image data
            [im1, units1] = imRead3D(self.im1Name);
            
            % Detect keypoints
            keys = detectSift3D(im1, 'units', units1);
            
            % Extract descriptors using the pyramid
            [descPyr, coordsPyr] = extractSift3D(keys);
            
            % Extract raw descriptors
            [descRaw, coordsRaw] = extractSift3D(keys, im1, units1);
            
            % Check the results
            assertElementsAlmostEqual(coordsPyr, coordsRaw);
            assertElementsAlmostEqual(descPyr, descRaw, 'absolute', 0.2);
            
        end
        
        % Test that "raw" keypoint orientations are close to those
        % extracted from a Gaussian scale-space pyramid
        function rawOrientationTest(self)
           if ~self.fullTest
              return 
           end
           
           % Load the image data
           [im, units] = imRead3D(self.im1Name);
           
           % Detect keypoints
           keys = detectSift3D(im, 'units', units);
           
           % Assign orientations to those same keypoints
           keysRaw = orientation3D(keys, im, units);
           
           % Check the dimensions
           assertEqual(size(keys), size(keysRaw));
           
           % Create a basis vector
           u = zeros(length(keys(1).ori), 1);
           u(1) = 1;
           
           % Check the results
           ang = zeros(length(keys), 1);
           for i = 1 : length(keys)
               key = keys(i);
               keyRaw = keysRaw(i);
               
               % Rotate the basis vector
               uRotKey = key.ori * u;
               uRotKeyRaw = keyRaw.ori * u;
               
               % Compute the angle between the rotated vectors
               ang(i) = acos(abs(dot(uRotKey, uRotKeyRaw)));
           end
           
           % Test the median angle
           assertElementsAlmostEqual(median(ang), 0, 'absolute', pi / 8);
        end
        
        % Test that detected keypoints are valid
        function detectValidTest(self)
            if ~self.fullTest
                return
            end
            
            % Load the image data
            [im, units] = imRead3D(self.im1Name);
            
            % Extract keypoints
            keys = detectSift3D(im, 'units', units);
            
            % Check the keypoints
            for i = 1 : length(keys)
                
                key = keys(i);
            
                % Check containment in the original image
                baseCoords = key.coords * pow2(-key.octave);
                assertTrue(all(baseCoords >= 0));
                assertTrue(all(baseCoords < size(im)));
                
                % Check orthogonality of the rotation matrix
                assertElementsAlmostEqual(key.ori * key.ori', ...
                    eye(length(key.ori)), 'absolute', 1E-3);
                
                % Check determinant of the rotation matrix
                assertElementsAlmostEqual(det(key.ori), 1, 'absolute', ...
                    1E-3);
            end
        end
        
        % Test registering an image against the CLI version
        function regCliTest(self)
            if ~self.fullTest || ispc
                return
            end
            
            % Output file names
            matchesName = 'matches.csv';
            transformName = 'transform.csv';
            
            % Register with the CLI
            status = runCmd([self.regCmd ' --matches ' matchesName ...
                ' --transform ' transformName ' ' self.im1Name ' ' ...
                self.im2Name]);
            assertEqual(status, 0);
            
            % Read the results
            matchesCli = csvread(matchesName);
            transformCli = csvread(transformName);
            
            % Convert the results to Matlab's format
            matchSrcCli = matchesCli(:, 1 : 3);
            matchRefCli = matchesCli(:, 4 : end);
            
            % Load the images
            [im1, units1] = imRead3D(self.im1Name);
            [im2, units2] = imRead3D(self.im2Name);
            
            % Register with Matlab
            [A, matchSrc, matchRef] = registerSift3D(im1, im2, ...
                'srcUnits', units1, 'refUnits', units2);
            
            % Check the dimensions
            assertEqual(size(A), size(transformCli));
            assertEqual(size(matchSrc), size(matchSrcCli));
            assertEqual(size(matchRef), size(matchRefCli));
            
            % Check the matches (the only error is conversion to text)
            assertElementsAlmostEqual(matchSrc, matchSrcCli, ...
                'absolute', self.tolText);
            assertElementsAlmostEqual(matchRef, matchRefCli, ...
                'absolute', self.tolText);
            
            % Check the transformation (discrepancies introduced by 
            % randomized regression)
            assertElementsAlmostEqual(A(:, 1 : 3), ...
                transformCli(:, 1 : 3), 'absolute', 5E-2);
            assertElementsAlmostEqual(A(:, end), transformCli(:, end), ...
                'absolute', 5);
            
            % Clean up
            delete(matchesName);
            delete(transformName);
        end
        
        % Test anisotropic registration
        function regAnisoTest(self)
            
            if ~self.fullTest
                return
            end
            
            % Load the image
            [im, units] = imRead3D(self.im1Name);
            
            % Remove half of the slices of the image
            imAniso = im(:, :, 1 : 2 : end);
            unitsAniso = [units(1) units(2) units(3) * 2];
            
            % Register the original to the anisotropic image
            A = registerSift3D(im, imAniso, 'srcUnits', units, ...
                'refUnits', unitsAniso, 'resample', true);
            
            % Form the reference (ground truth) transformation
            refA = [eye(3) zeros(3, 1)];
            refA(3, 3) = 2;
            
            % Check the transformation
            assertElementsAlmostEqual(A(:, 1 : 3), refA(:, 1 : 3), ...
                'absolute', 5E-2);
            assertElementsAlmostEqual(A(:, end), refA(:, end), ...
                'absolute', 5);
        end
        
        % Test matching descriptors against the C version
        function matchTest(self)
            
            if ~self.fullTest
                return
            end
            
            % Load the images
            [im1, units1] = imRead3D(self.im1Name);
            [im2, units2] = imRead3D(self.im2Name);
            
            % Register with the C version and get the matches
            [~, match1, match2] = registerSift3D(im1, im2, ...
                'srcUnits', units1, 'refUnits', units2);
            
            % Extract descriptors from each image
            keys = detectSift3D(im1, 'units', units1);
            [desc1, coords1] = extractSift3D(keys);
            keys = detectSift3D(im2, 'units', units2);
            [desc2, coords2] = extractSift3D(keys);
            
            % Match with the Matlab version and convert to coordinates
            matches = matchSift3D(desc1, coords1, desc2, coords2);
            match1M = coords1(matches(:, 1), :);
            match2M = coords2(matches(:, 2), :);
            
            assertElementsAlmostEqual(match1M, match1, 'relative', 1E-3);
            assertElementsAlmostEqual(match2M, match2, 'relative', 1E-3);
            
        end
        
        % Test registration with invalid matching threshold
        function regInvalidMatchTest(self)
            
           % Load the images
           im1 = imRead3D(self.im1Name);
           im2 = imRead3D(self.im2Name);
            
           threwErr = false;
           try 
               A = registerSift3D(im1, im2, 'nnThresh', 2);
           catch ME
               threwErr = true;
           end
           assertTrue(threwErr);
        end
        
        % Test registration with invalid error threshold
        function regInvalidErrTest(self)
            
            % Load the images
            im1 = imRead3D(self.im1Name);
            im2 = imRead3D(self.im2Name);
            
            threwErr = false;
            try
                A = registerSift3D(im1, im2, 'errThresh', -1);
            catch ME
                threwErr = true;
            end
            assertTrue(threwErr);
        end
        
        % Test registration with invalid number of iterations
        function regInvalidIterTest(self)
            
            % Load the images
            im1 = imRead3D(self.im1Name);
            im2 = imRead3D(self.im2Name);
            
            threwErr = false;
            try
                A = registerSift3D(im1, im2, 'numIter', 0);
            catch ME
                threwErr = true;
            end
            assertTrue(threwErr);
        end
        
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
            
            % Remove any past instances of this file
            if exist(imName, 'file')
                delete(imName)
            end
            
            % Make random image data, scaled to the range [0, 1]
            imWritten = rand(10, 15, 20);
            imWritten = imWritten / max(imWritten(:));
            
            % Write the image as a DICOM file
            imWrite3D(imName, imWritten);
            
            % Read the image back and scale it
            imRead = imRead3D(imName);
            imRead = imRead / max(imRead(:));
            
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
            
            % Read the image back and scale it
            imRead = imRead3D(dirName);
            imRead = imRead / max(imRead(:));
            
            % Clean up
            rmdir(dirName, 's');
            
            % Ensure the results are identical
            assertElementsAlmostEqual(imWritten, imRead, 'absolute', 1E-2);
        end
        
        % Test reading and writing a 2D NIFTI image
        function nifti2DTest(self)
            % The temporary file name
            imName = 'temp.nii.gz';
            
            % Make random image data
            imWritten = rand(20, 15);
            
            % Write the image as a NIFTI file
            imWrite3D(imName, imWritten);
            
            % Read the image back
            imRead = imRead3D(imName);
            
            % Clean up
            delete(imName);
            
            % Ensure the results are identical
            assertElementsAlmostEqual(imWritten, imRead, 'relative', 1E-3);
        end
        
        % Test reading and writing a 2D DICOM image
        function dicom2DTest(self)
            
            % The temporary file name
            imName = 'temp.dcm';
            
            % Make random image data, scaled to the range [0, 1]
            imWritten = rand(20, 15);
            imWritten = imWritten / max(imWritten(:));
            
            % Write the image as a DICOM file
            imWrite3D(imName, imWritten);
            
            % Read the image back and scale it
            imRead = imRead3D(imName);
            imRead = imRead / max(imRead(:));
            
            % Clean up
            delete(imName);
            
            % Ensure the results are identical
            assertElementsAlmostEqual(imWritten, imRead, 'absolute', 1E-2);
        end
        
        % Test reading and writing units from a NIFTI image
        function niftiUnitsTest(self)
            
            % The temporary file name
            imName = 'temp.nii.gz';
            
            % Make random image data
            imWritten = rand(10, 15, 20);
            
            % Make random units
            unitsWritten = rand(3, 1);
            
            % Write the image as a NIFTI file
            imWrite3D(imName, imWritten, unitsWritten);
            
            % Read the units back
            [~, unitsRead] = imRead3D(imName);
            
            % Clean up
            delete(imName);
            
            % Ensure the results are identical
            assertElementsAlmostEqual(unitsWritten, unitsRead, ...
                'relative', 1E-3);
        end
        
        % Test reading and writing units from a DICOM image
        function dicomUnitsTest(self)
            
            % The temporary file name
            imName = 'temp.dcm';
            
            % Make random image data
            imWritten = rand(10, 15, 20);
            
            % Make random units
            unitsWritten = rand(3, 1);
            
            % Write the image as a NIFTI file
            imWrite3D(imName, imWritten, unitsWritten);
            
            % Read the units back
            [~, unitsRead] = imRead3D(imName);
            
            % Clean up
            delete(imName);
            
            % Ensure the results are identical
            assertElementsAlmostEqual(unitsWritten, unitsRead, ...
                'relative', 1E-3);
        end
        
        % Test reading and writing units from a directory of DICOM images
        function dirUnitsTest(self)
            
            % The temporary file name
            imName = 'temp';
            
            % Make random image data
            imWritten = rand(10, 15, 20);
            
            % Make random units
            unitsWritten = rand(3, 1);
            
            % Write the image as a NIFTI file
            imWrite3D(imName, imWritten, unitsWritten);
            
            % Read the units back
            [~, unitsRead] = imRead3D(imName);
            
            % Clean up
            rmdir(imName, 's');
            
            % Ensure the results are identical
            assertElementsAlmostEqual(unitsWritten, unitsRead, ...
                'relative', 1E-3);
        end
        
        % Test reading and writing units from a 2D DICOM image
        function units2DTest(self)
            
            % The temporary file name
            imName = 'temp.dcm';
            
            % Make random image data
            imWritten = rand(20, 15);
            
            % Make random units
            unitsWritten = rand(2, 1);
            
            % Write the image as a NIFTI file
            imWrite3D(imName, imWritten, unitsWritten);
            
            % Read the image back
            [~, unitsRead] = imRead3D(imName);
            
            % Remove the trailing units
            unitsRead = unitsRead(1 : length(unitsWritten));
            
            % Clean up
            delete(imName);
            
            % Ensure the results are identical
            assertElementsAlmostEqual(unitsWritten, unitsRead, ...
                'relative', 1E-3);
        end
        
        % Test switching the parameter order in imWrite3D
        function writeSwappedParamsTest(self)
            
            % Make a random image
            im = rand(10, 10, 10);
            
            % Try to write it, with parameters exchanged
            threwErr = false;
            try
                imWrite3D(im, 'im.nii.gz');
            catch ME
                threwErr = true;
            end
            assertTrue(threwErr);
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
        
        % Test writing negative units
        function negativeUnitsTest(self)
            
            % Invalid units
            units = [1 1 -1];
            
            % Fake image
            im = zeros(10, 10, 10);
            
            threwErr = false;
            try
                imWrite3D('fake.nii.gz', im, units);
            catch ME
                threwErr = true;
            end
            assertTrue(threwErr);
        end
        
        % Test writing zero-valued units
        function zeroUnitsTest(self)
            
            % Invalid units
            units = [1 0 1];
            
            % Fake image
            im = zeros(10, 10, 10);
            
            threwErr = false;
            try
                imWrite3D('fake.nii.gz', im, units);
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

% The CLI is not supported on Windows
if ispc
    warning(['The command-line interface is not supported in the ' ...
        'Windows version of SIFT3D'])
end

% Strip the LD_LIBRARY_PATH environment variable of Matlab
% directories
ldPathVar = 'LD_LIBRARY_PATH';
oldLdPath = getenv(ldPathVar);
newLdPath = regexprep(oldLdPath, '[^:]*MATLAB[^:]*:*', '', 'ignorecase');
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
