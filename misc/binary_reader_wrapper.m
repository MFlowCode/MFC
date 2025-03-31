% Wrapper for reading MFC binary output file
% Works for 1D/2D/3D with multiple processors

function pres = binary_reader_wrapper(binDir, ti, tf, t_delta, dim)
    % Add more output variables (like rho or xCoords) as desired

    % Total time steps
    tArr = ti : t_delta : tf;
    tArrLen = length(tArr);

    if (dim ~= 1)
        [iProcList, m, n, p, xIdxs, yIdxs, zIdxs, xCoords, yCoords, zCoords] = getProcIdx(binDir, tArr, dim);
    else
        % For 1D, the folder 'root' contains all the data
        if (isempty(binDir))
            error(strcat("ERROR: invalid binDir (not a binary/root folder): ", binDir))
        end
        filename = fullfile(binDir, 'root', [num2str(tArr(1)), '.dat']);
        dat = f_binary_reader(filename, 'n', 'real*8', 50);
        m = dat.m;
        iProcList = 1;
        xIdxs{1} = 1:m+1;
    end

    % Loop through files for each time step
    for tIdx = 1:tArrLen
        if (mod(tIdx, 10) == 0 || tIdx == 1)
            disp(['Reading time step ', num2str(tIdx), ' of ', num2str(tArrLen)]);
        end
        
        for iProc = 1 : length(iProcList)
            if (dim ~= 1)
                filename = fullfile(binDir, ['p', num2str(iProcList(iProc)-1)], [num2str(tArr(tIdx)), '.dat']);
            else
                filename = fullfile(binDir, 'root', [num2str(tArr(tIdx)), '.dat']);
            end
            dat = f_binary_reader(filename, 'n', 'real*8', 50);

            xIdx = xIdxs{iProc};
            nx = m(iProc)+1;
            if (dim >= 2)
                yIdx = yIdxs{iProc};
                ny = n(iProc)+1;
            end
            if (dim == 3)
                zIdx = zIdxs{iProc};
                nz = p(iProc)+1;
            end
            
            % Add more variables (like rho) as desired
            if (dim == 1)
                pres(xIdx, tIdx) = dat.pres;
            elseif (dim == 2)
                pres(xIdx, yIdx, tIdx) = reshape(dat.pres, nx, ny);
            elseif (dim == 3)
                pres(xIdx, yIdx, zIdx, tIdx) = reshape(dat.pres, nx, ny, nz);
            end
        end
    end
end

%% Helper Functions

function [iProcList, m, n, p, xIdxs, yIdxs, zIdxs, xCoords, yCoords, zCoords] = getProcIdx(binDir, tArr, dim)
% Set up global index mapping for multiple processors (only for 2D/3D, so dim == 2 or 3)

% Returns:
%   iProcList - list of folder indices corresponding to valid processors
%   m, n, p   - number of cells in each proc (lists)
%   xIdxs, yIdxs, zIdxs       - global index arrays for each proc
%   xCoords, yCoords, zCoords - global cell center coordinate arrays

    % List folders according to processor number used for the simulation
    p_folders = dir( fullfile(binDir, 'p*') ); 
    nProcFolders = length(p_folders);
    if (nProcFolders == 0)
        disp("ERROR: No p_* folders found in binDir:")
        disp(binDir)
        error('No processor folders found.');
    end

    % First get the coord range for each proc (using the first time step)
    iProcList = [];
    validProc = 0;
    for procNum = 1:nProcFolders
        filename = fullfile(binDir, ['p', num2str(procNum-1)], [num2str(tArr(1)), '.dat']);
        dat = f_binary_reader(filename, 'n', 'real*8', 50);
        if ((dat.m == 0) || (dim >= 2 && dat.n == 0) || (dim == 3 && dat.p == 0))
            continue
        end
        validProc = validProc + 1;
        m(validProc) = dat.m;
        n(validProc) = dat.n;
        xCoord{validProc} = dat.x_cb;
        yCoord{validProc} = dat.y_cb;
        if (dim == 3)
            p(validProc) = dat.p;
            zCoord{validProc} = dat.z_cb;
        end
        iProcList(end+1) = procNum;
    end
    nProc = validProc;
    assert(nProc == length(iProcList));

    % Assign empty arrays if not used
    if dim < 3
        p = [];
        zCoord = {};
    end

    % For each proc we build a list of processors that come before it 
        % in a given direction by comparing the coordinates
    xProcOrder = cell(1, nProc);
    yProcOrder = cell(1, nProc);
    zProcOrder = {};
    if (dim == 3)
        zProcOrder = cell(1, nProc);
    end

    for iProc = 1 : nProc
        if (dim == 2)
            xProcOrder{iProc} = getProcOrder(xCoord, iProc, yCoord);
            yProcOrder{iProc} = getProcOrder(yCoord, iProc, xCoord);
        elseif (dim == 3)
            xProcOrder{iProc} = getProcOrder(xCoord, iProc, yCoord, zCoord);
            yProcOrder{iProc} = getProcOrder(yCoord, iProc, xCoord, zCoord);
            zProcOrder{iProc} = getProcOrder(zCoord, iProc, xCoord, yCoord);
        end
    end

    % Using the order, we compute the global indices using the proc orders
    for iProc = 1 : nProc
        xLocal{iProc} = 1:(m(iProc)+1);
        yLocal{iProc} = 1:(n(iProc)+1);
        if (dim == 3)
            zLocal{iProc} = 1:(p(iProc)+1);
        end
    end

    xIdxs = getGlobalIdx(xLocal, xProcOrder);
    yIdxs = getGlobalIdx(yLocal, yProcOrder);
    zIdxs = {};
    if (dim == 3)
        zIdxs = getGlobalIdx(zLocal, zProcOrder);
    end

    % Get global cell center locations
    for iProc = 1 : nProc
        xCoords( xIdxs{iProc} ) = xCoord{iProc}(1:end-1) + diff(xCoord{iProc})/2;
        yCoords( yIdxs{iProc} ) = yCoord{iProc}(1:end-1) + diff(yCoord{iProc})/2;
        if (dim == 3)
            zCoords( zIdxs{iProc} ) = zCoord{iProc}(1:end-1) + diff(zCoord{iProc})/2;
        end
    end
    % Assign empty outputs if not used
    if dim < 3
        zCoords = [];
    end
end

function order = getProcOrder(coordCell, iProc, varargin)
% For the processor, compute a list of processors that come it based on their coordinates

    order = [];
    currentCoord = coordCell{iProc}(1);
    nProc = numel(coordCell);
    for j = 1:nProc
        if j == iProc
            continue;
        end
        if coordCell{j}(1) < currentCoord
            valid = true;
            for k = 1:length(varargin)
                if varargin{k}{iProc}(1) ~= varargin{k}{j}(1)
                    valid = false;
                    break;
                end
            end
            if valid
                order(end+1) = j;
            end
        end
    end
end

function globalIdxs = getGlobalIdx(localIdxs, procOrder)
% Compute the global indices for each processor

    nProc = length(localIdxs);
    globalIdxs = cell(1, nProc);
    for i = 1:nProc
        offset = 0;
        if ~isempty(procOrder{i})
            for j = procOrder{i}
                offset = offset + length(localIdxs{j});
            end
        end
        globalIdxs{i} = localIdxs{i} + offset;
    end
end