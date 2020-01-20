%% get_unique_ids

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% October 2019; Last revision:

function unique_id = get_unique_ids(animal_ids)

unique_id = cumsum([1 diff([animal_ids])~=0]);
