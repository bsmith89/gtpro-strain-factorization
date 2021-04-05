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

CREATE TABLE mgen
  ( mgen_id PRIMARY KEY
  , sample_id REFERENCES sample(sample_id)
  , plate_number
  , plate_well
  , filename_r1
  , filename_r2
  , mgen_notes
  );

CREATE TABLE mgen_x_mgen_group
  ( mgen_id REFERENCES mgen(mgen_id)
  , mgen_group
  );
