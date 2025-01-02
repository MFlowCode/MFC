% @brief The following function takes care of reading the binary database
%			   files that may be outputted by the MFC during the post-processing
%			   of the raw data files. The function receives the name of a data-
%			   base file, its location, byte order and precision it was written
%			   with, and finally, the maximum allowable length of the database's
%			   variables' names. The function then reads the data stored in the
%			   database file and generates an equivalent MATLAB structure array.
%			   Note that the byte order and the precision inputs are in the same
%			   notation and feature the same options as those offered by native
%			   MATLAB functions.
function mat_struct = f_binary_reader(file_loc, byte_order, precision, name_len)
	
	% NOTE: The user may notice that throughout this file the following command
	% is frequently invoked - fread(dbfile,1,'real*4',byte_order). This is made
	% necessary by the unformatted Fortran output that is used to generate the
	% binary database files. The purpose of this command is to read and discard
	% the purposeless header and footer information that divides the individual
	% records located in the database.
	
	
	% Opening the binary database file and generating its handle
	dbfile = fopen(file_loc, 'r', byte_order);
	
	
	% Input and Output of the Binary Database Structure Information
	
	% Reading in, respectively, cell-numbers in the x-, y- and z-directions and
	% the total number of flow variables located in the database
	fread(dbfile, 1, 'real*4', byte_order);
	m      = fread(dbfile, 1, 'integer*4', byte_order);
	n      = fread(dbfile, 1, 'integer*4', byte_order);
	p      = fread(dbfile, 1, 'integer*4', byte_order);
	dbvars = fread(dbfile, 1, 'integer*4', byte_order);
	fread(dbfile, 1, 'real*4', byte_order);
	
	% Writing the cell-numbers in the x-, y- and z-directions to output
	mat_struct.('m') = m;
	mat_struct.('n') = n;
	mat_struct.('p') = p;
	
	% Input and Output of the Grid Data
	
	% Reading in the cell-boundary locations in the x-direction
	fread(dbfile, 1, 'real*4', byte_order);
	x_cb = fread(dbfile, m+2, precision, byte_order);
	
	% Writing the cell-boundary locations in the x-direction to output
	mat_struct.('x_cb') = x_cb;
	
	if(n > 0)
		
		% Reading in the cell-boundary locations in the y-direction
		y_cb = fread(dbfile, n+2, precision, byte_order);
		
		% Writing the cell-boundary locations in the y-direction to output
		mat_struct.('y_cb') = y_cb;
		
		if(p > 0)
			
			% Reading in the cell-boundary locations in the z-direction
			z_cb = fread(dbfile, p+2, precision, byte_order);
			
			% Writing cell-boundary locations in the z-direction to output
			mat_struct.('z_cb') = z_cb;
			
			% Setting size of the grid for a 3D simulation
			grid_size = [m+1, n+1, p+1];
			
		else
			
			% Setting size of the grid for a 2D simulation
			grid_size = [m+1, n+1];
			
		end
		
	else
		
		% Setting size of the grid for a 1D simulation
		grid_size = m+1;
		
	end
	
	% Discarding the footer of the grid data record
	fread(dbfile, 1, 'real*4', byte_order);
	
	% Input and Output of the Flow Variables
	for j = 1:dbvars
		
		% Reading in a flow variable from the database
		fread(dbfile, 1, 'real*4', byte_order);
		varname = strtrim(fread(dbfile, name_len, '*char', byte_order)');
		q_fp = fread(dbfile, prod(grid_size), precision, byte_order);
		fread(dbfile, 1, 'real*4', byte_order);
		
		% Writing the flow variable to output
		mat_struct.(varname) = q_fp;
		
	end
	
	% Closing the binary database file
	fclose(dbfile);
	
% END: function f_mfc_binary_reader
