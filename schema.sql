CREATE TABLE subject
  ( subject_id PRIMARY KEY
  , subject_notes
  );

CREATE TABLE sample
  ( sample_id PRIMARY KEY
  , subject_id REFERENCES subject(subject_id)
  , collection_date DATE
  , storage_location
  , sample_weight
  , sample_notes
  );

CREATE TABLE mgen_library
  ( mgen_library_id PRIMARY KEY
  , sample_id REFERENCES sample(sample_id)
  , plate_number
  , plate_well
  , filename_r1
  , filename_r2
  , mgen_library_notes
  );

CREATE TABLE mgen_library_x_mgen_library_group
  ( mgen_library_id REFERENCES mgen_library(mgen_library_id)
  , mgen_library_group
  );
