-- Stop at the first error, such as a missing CSV file
\set ON_ERROR_STOP on

DROP TABLE IF EXISTS activsg2000;
CREATE TABLE activsg2000 (
    sample   bigint,
    t        double precision,
    class    text,
    id       text,
    variable text,
    value    double precision,
    -- t repeats at events, so sample identifies the time point
    PRIMARY KEY (sample, class, id, variable)
);
\copy activsg2000 FROM 'ACTIVSg2000.long.csv' WITH (FORMAT csv, HEADER)

-- Generators with the largest speed deviation
SELECT class, id, ROUND(MAX(ABS(value))::numeric, 6) AS max_speed_deviation
FROM activsg2000
WHERE variable = 'omega'
GROUP BY class, id
ORDER BY max_speed_deviation DESC
LIMIT 5;
