SELECT a.patientTableId, a.patientName, b.patientTableId, b.studyId, b.studyFileDir, b.studyDate
FROM patient a
LEFT JOIN study b
on a.patientTableId = b.patientTableId
WHERE b.studyDate > '2018-01-01'

WITH study_result AS(
SELECT a.patientTableId, a.patientName, b.patientTableId, b.studyId, b.studyFileDir, b.studyDate
FROM patient a
LEFT JOIN study b
on a.patientTableId = b.patientTableId
WHERE b.studyDate > '2017-01-26' AND a.patientName = 'KEPLER^JUDITH'
),
series_result AS(
SELECT * FROM series_all
WHERE studyTableId IN (SELECT studyId FROM study_result)
)

SELECT * FROM series_result

WITH study_result AS(
SELECT a.patientTableId, a.patientName, b.patientTableId, b.studyId, b.studyFileDir, b.studyDate
FROM patient a
LEFT JOIN study b
on a.patientTableId = b.patientTableId
WHERE b.studyDate > '2017-01-26' AND a.patientName = 'KEPLER^JUDITH'
)


SELECT a.studyId, a.patientName, a.studyDate, b.kspaceDataTableId, b.kspaceTableID, b.kspaceFileDir, b.studyTableId
FROM study_result a
LEFT JOIN kspace b
on a.studyId = b.studyTableId


conn = sqlite('./Copy_of_db/test_db.sqlite3', 'readonly');
sqlquery = ...
    ['select distinct reconstruction.kspaceId, reconstruction.kspaceId2, reconstruction.referenceFileDir, motion_signal.signal, reconstruction.timestep ' ...
    'from reconstruction, states, dvfs, motion_signal ' ...
    'where reconstruction.reconstructionIndex = states.reconstructionIndex and states.reconstructionIndex = 2 and states.signalType = "b" and dvfs.reconstructionIndex = reconstruction.reconstructionIndex ' ...
    'and motion_signal.reconstructionIndex = reconstruction.reconstructionIndex and motion_signal.signalType = "b"'];
disp(sqlquery);
results = fetch(conn,sqlquery);
disp(results);
close(conn);

kSpaceID = [results.kspaceId, results.kspaceId2]; % in db
referenceFileDir = results.referenceFileDir;
scanDir = convertStringsToChars(append("/mnt/ibrixfs01-FUNCI/", extractBetween(referenceFileDir,"/nfs/corexfs/FUNCI/", "/reconstruction")));
reconstructionDir = fullfile(scanDir, 'reconstruction');
zsi = results.signal.split(',');
zSI = str2double(zsi);

conn = sqlite('./Copy_of_db/test_db.sqlite3', 'readonly');
sqlquery4 = ...
    ['select distinct reconstruction.kspaceId, reconstruction.kspaceId2, reconstruction.referenceFileDir, motion_signal.signal, reconstruction.timestep ' ...
    'from reconstruction, states, dvfs, motion_signal ' ...
    'where reconstruction.reconstructionIndex = states.reconstructionIndex and states.reconstructionIndex = 2 and states.signalType = "b" and dvfs.reconstructionIndex = reconstruction.reconstructionIndex ' ...
    'and motion_signal.reconstructionIndex = reconstruction.reconstructionIndex and motion_signal.signalType = "g"'];
disp(sqlquery4);
results4 = fetch(conn,sqlquery4);
disp(results4);
close(conn);

phi = results4.signal.split(',');
phi = str2double(phi);

conn = sqlite('./Copy_of_db/test_db.sqlite3', 'readonly');
sqlquery1 = ...
['select dvfs.defIndex, dvfs.fileDir '...
'from reconstruction, dvfs '...
'where reconstruction.reconstructionIndex = dvfs.reconstructionIndex and reconstruction.reconstructionIndex = 2 and dvfs.signalType = "b"'];
disp(sqlquery1);
results1 = fetch(conn,sqlquery1);
disp(results1);
close(conn);

b_state_num = results1.defIndex(end);
b_dvfs = results1.fileDir;

currentTransforms = cell(1, b_state_num);

for iMotionState = 1:b_state_num
    currentTransforms{iMotionState} = loadArrayVolumeFromNii(convertStringsToChars(append("/mnt/ibrixfs01-FUNCI/", extractAfter(b_dvfs(iMotionState),"/nfs/corexfs/FUNCI/")))); % the input of this function is just a file dir, so can be directly load from the db
end

conn = sqlite('./Copy_of_db/test_db.sqlite3', 'readonly');
sqlquery2 = ...
['select dvfs.defIndex, dvfs.fileDir '...
'from reconstruction, dvfs '...
'where reconstruction.reconstructionIndex = dvfs.reconstructionIndex and reconstruction.reconstructionIndex = 2 and dvfs.signalType = "s"'];
disp(sqlquery2);
results2 = fetch(conn,sqlquery2);
disp(results2);
close(conn);

s_state_num = results2.defIndex(end);
s_dvfs = results2.fileDir;

currentTransforms2 = cell(1, s_state_num);

for iMotionState = 1:s_state_num
    currentTransforms2{iMotionState} = loadArrayVolumeFromNii(convertStringsToChars(append("/mnt/ibrixfs01-FUNCI/", extractAfter(s_dvfs(iMotionState),"/nfs/corexfs/FUNCI/")))); 
end

conn = sqlite('./Copy_of_db/test_db.sqlite3', 'readonly');
sqlquery3 = ...
['select dvfs.defIndex, dvfs.fileDir '...
'from reconstruction, dvfs '...
'where reconstruction.reconstructionIndex = dvfs.reconstructionIndex and reconstruction.reconstructionIndex = 2 and dvfs.signalType = "g"'];
disp(sqlquery3);
results3 = fetch(conn,sqlquery3);
disp(results3);
close(conn);

g_state_num = results3.defIndex(end);
g_dvfs = results3.fileDir;

currentTransforms3 = cell(1, g_state_num);

for iMotionState = 1:g_state_num
    currentTransforms3{iMotionState} = loadArrayVolumeFromNii(convertStringsToChars(append("/mnt/ibrixfs01-FUNCI/", extractAfter(g_dvfs(iMotionState),"/nfs/corexfs/FUNCI/")))); 
end

referenceState = loadArrayVolumeFromNii(convertStringsToChars(append("/mnt/ibrixfs01-FUNCI/", extractAfter(referenceFileDir,"/nfs/corexfs/FUNCI/"))));

timestep = results.timestep;


% suppose from 4000-4408 as a 60s period
startspokes = 4000;

%breathing
% take the exhale state
[~, iMotionStates1] = sort(zSI);    
iMotionStates1(iMotionStates1) = 1:numel(iMotionStates1);
% iState1 = iMotionStates1(N)*21/7000; % get non-integer state number
iState1 = b_state_num;
transform_b = currentTransforms{iState1};
transform_b_inv = invertDeformation(referenceState, transform_b);

%slow configuration
transform_s_invs = cell(1, 60);
t = 1;
t1 = 1;
spoke_index = 0;
while (t < 61)
    if (t1 >= t) 
        iState2 = (startspokes + spoke_index) / 7000 * double(s_state_num); % total number of spoke should be put into db
        temp1 = reshape(currentTransforms2{floor_in_house(iState2)+1}.A, 1, []);
        temp2 = reshape(currentTransforms2{ceil_in_house(iState2)+1}.A, 1, []);
        interped = interp1([floor_in_house(iState2)+1  ceil_in_house(iState2)+1], [temp1', temp2']', iState2+1, 'linear');
        transform_s = cloneWithNewContent(currentTransforms2{1}, reshape(interped, size(currentTransforms2{1}.A)));
        transform_s_inv = invertDeformation(referenceState, transform_s);
        transform_s_invs{t} = transform_s_inv; % 1 - 60
        outputfileName = sprintf('output3/slow_def_%d.img', t);
        exportArrayVolumeToNii(transform_s_inv, outputfileName);
        finalVolume = resampleToDeformationField(referenceState, transform_s_inv, 'spline', 'none'); %final reconstruction by deforming reference image
        outputfileName = sprintf('output3/slow_state_%d.img', t);
        exportArrayVolumeToNii(finalVolume, outputfileName);
        t = t + 1;
    end
    t1 = t1 + timestep;
    spoke_index = spoke_index + 1;
end


%gi motion
[~, iMotionStates3] = sort(phi);    
iMotionStates3(iMotionStates3) = 1:numel(iMotionStates3);
% transform_g_invs = cell(1, 60);
t = 1;
t1 = 1;
spoke_index = 0;
while (t < 61)
    if (t1 >= t) 
        iState3 = iMotionStates3(startspokes + spoke_index) / 7000 * double(g_state_num);
        if iState3 > double(g_state_num) - 1
            iState3 = double(g_state_num) - 0.001;
        end
        fprintf('%d', iState3);
        temp1 = reshape(currentTransforms3{floor_in_house(iState3)+1}.A, 1, []);
        temp2 = reshape(currentTransforms3{ceil_in_house(iState3)+1}.A, 1, []);
        interped = interp1([floor_in_house(iState3)+1  ceil_in_house(iState3)+1], [temp1', temp2']', iState3+1, 'linear');
        transform_g = cloneWithNewContent(currentTransforms3{1}, reshape(interped, size(currentTransforms3{1}.A)));
        transform_g_inv = invertDeformation(referenceState, transform_g);
        defFieldArray_1 = resampleToDeformationField(transform_s_invs{t},transform_g_inv,'linear');
        defFieldArray = resampleToDeformationField(defFieldArray_1,transform_b_inv,'linear'); %combine 3 deformations
        finalVolume = resampleToDeformationField(referenceState, defFieldArray, 'spline', 'none'); %final reconstruction by deforming reference image
        outputfileName = sprintf('output3/gi_state_%d_%f.img', t, phi(startspokes + spoke_index));
        exportArrayVolumeToNii(finalVolume, outputfileName);
        outputfileName = sprintf('output3/gi_def_%d.img', t);
        exportArrayVolumeToNii(defFieldArray, outputfileName); % this is not the gi def, but the def from final reference state to gi
        t = t + 1;
    end
    t1 = t1 + timestep;
    spoke_index = spoke_index + 1;
end




function result = floor_in_house(num)
if ceil(num) == floor(num)
    result = num;
else
    result = floor(num);
end
end

function result = ceil_in_house(num)
if ceil(num) == floor(num)
    result = num + 1;
else
    result = ceil(num);
end
end