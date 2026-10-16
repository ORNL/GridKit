\set ON_ERROR_STOP on
BEGIN;
DROP TABLE IF EXISTS activsg2000;
CREATE TABLE activsg2000 (
    sample   bigint           NOT NULL,
    t        double precision NOT NULL,
    class    text             NOT NULL,
    id       text             NOT NULL,
    variable text             NOT NULL,
    value    double precision NOT NULL,
    PRIMARY KEY (sample, class, id, variable)
);
\copy activsg2000 FROM 'ACTIVSg2000.long.csv' WITH (FORMAT csv, HEADER)
COMMIT;

SELECT
    class,
    id,
    ROUND(MAX(ABS(value))::numeric, 6) AS max_speed_deviation
FROM activsg2000
WHERE variable = 'omega'
GROUP BY class, id
ORDER BY max_speed_deviation DESC
LIMIT 5;
