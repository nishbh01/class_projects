-- Removing tables with same name as tables being created
DROP TABLE DOCTOR CASCADE CONSTRAINTS;
DROP TABLE PATIENT CASCADE CONSTRAINTS;
DROP TABLE PRESCRIPTION CASCADE CONSTRAINTS;
DROP TABLE SIDE_EFFECTS CASCADE CONSTRAINTS;

-- Creating tables for Health Care Datbase
CREATE TABLE doctors (
id INTEGER PRIMARY KEY,
first_name VARCHAR(50) NOT NULL,
last_name VARCHAR(50) NOT NULL,
email VARCHAR(100) NOT NULL UNIQUE,
phone_number VARCHAR(20) NOT NULL,
specialization VARCHAR(100) NOT NULL,
address_line1 VARCHAR(100) NOT NULL,
address_line2 VARCHAR(100),
city VARCHAR(50) NOT NULL,
state VARCHAR(50) NOT NULL,
zip_code VARCHAR(10) NOT NULL
);


CREATE TABLE patients (
id INTEGER PRIMARY KEY,
first_name VARCHAR(50) NOT NULL,
last_name VARCHAR(50) NOT NULL,
email VARCHAR(100) NOT NULL,
phone_number VARCHAR(20) NOT NULL,
gender VARCHAR(10) NOT NULL,
birthdate DATE NOT NULL,
address_line1 VARCHAR(100) NOT NULL,
address_line2 VARCHAR(100),
city VARCHAR(50) NOT NULL,
state VARCHAR(50) NOT NULL,
zip_code VARCHAR(10) NOT NULL,
doc_id INTEGER,
FOREIGN KEY (doc_id) REFERENCES doctors(id)
);


CREATE TABLE appointments (
id INTEGER PRIMARY KEY,
doctor_id INTEGER NOT NULL,
patient_id INTEGER NOT NULL,
date_time TIMESTAMP NOT NULL,
duration INTEGER NOT NULL,
reason VARCHAR(200),
status VARCHAR(20) NOT NULL CHECK (status IN ('scheduled', 'completed', 'cancelled')),
FOREIGN KEY (doctor_id) REFERENCES doctors(id),
FOREIGN KEY (patient_id) REFERENCES patients(id)
);


CREATE TABLE medications (
id INTEGER PRIMARY KEY,
name VARCHAR(100) NOT NULL,
price DECIMAL(10, 2) NOT NULL,
side_effects VARCHAR(500),
instructions VARCHAR(500),
dosage VARCHAR(100),
category VARCHAR(50),
pat_id INTEGER,
FOREIGN KEY (pat_id) REFERENCES patients(id)
);

-- Inserting values in each entity
INSERT INTO doctors VALUES (1, 'John', 'Smith', 'johnsmith@gmail.com', '555-123-4567', 'Cardiologist', '123 Main St', NULL, 'Los Angeles', 'CA', '90001');
INSERT INTO doctors VALUES(2, 'Jane', 'Doe', 'janedoe@gmail.com', '555-987-6543', 'Dermatologist', '456 Elm St', NULL, 'New York', 'NY', '10001');
INSERT INTO doctors VALUES(3, 'Bob', 'Johnson', 'bjohnson@yahoo.com', '555-555-1212', 'Pediatrician', '789 Oak Ave', NULL, 'Chicago', 'IL', '60601');
INSERT INTO doctors VALUES(4, 'Susan', 'Lee', 'susanlee@gmail.com', '555-555-2323', 'Neurologist', '321 Pine St', NULL, 'San Francisco', 'CA', '94101');
INSERT INTO doctors VALUES(5, 'David', 'Kim', 'dkim@hotmail.com', '555-555-3434', 'Oncologist', '2345 Maple Ave', NULL, 'Houston', 'TX', '77001');


INSERT INTO patients VALUES(1, 'Jennifer', 'Brown', 'jenniferbrown@gmail.com', '555-555-1212', 'female', TO_DATE('1990-01-01', 'YYYY-MM-DD'), '123 Main St', NULL, 'Los Angeles', 'CA', '90001', 2);
INSERT INTO patients VALUES(2, 'Mike', 'Johnson', 'mjohnson@yahoo.com', '555-555-2323', 'male', TO_DATE('1985-02-02', 'YYYY-MM-DD'), '456 Elm St', NULL, 'New York', 'NY', '10001', 1);
INSERT INTO patients VALUES(3, 'Lisa', 'Lee', 'lisalee@gmail.com', '555-555-3434', 'female', TO_DATE('1975-03-03', 'YYYY-MM-DD'), '789 Oak Ave', NULL, 'Chicago', 'IL', '60601', 3);
INSERT INTO patients VALUES(4, 'Tom', 'Smith', 'tomsmith@hotmail.com', '555-555-4545', 'male', TO_DATE('1965-04-04', 'YYYY-MM-DD'), '321 Pine St', NULL, 'San Francisco', 'CA', '94101', 4);
INSERT INTO patients VALUES(5, 'Emily', 'Kim', 'emilykim@gmail.com', '555-555-5656', 'female', TO_DATE('1955-05-05', 'YYYY-MM-DD'), '2345 Maple Ave', NULL, 'Houston', 'TX', '77001', 5);


INSERT INTO appointments VALUES(1, 1, 1, TIMESTAMP '2023-05-01 10:00:00', 60, 'Chest pain', 'Aspirin', 'completed');
INSERT INTO appointments VALUES(2, 2, 2, TIMESTAMP '2023-05-02 11:30:00', 30, 'Acne', 'Tretinoin', 'completed');
INSERT INTO appointments VALUES(3, 3, 3, TIMESTAMP '2023-05-03 13:00:00', 45, 'Flu', 'Tamiflu', 'scheduled');
INSERT INTO appointments VALUES(4, 4, 4, TIMESTAMP '2023-05-04 14:30:00', 60, 'Headache', 'Ibuprofen', 'scheduled');
INSERT INTO appointments VALUES(5, 5, 5, TIMESTAMP '2023-05-05 16:00:00', 30, 'Cancer', 'Chemotherapy', 'scheduled');


INSERT INTO medications VALUES(1, 'Lisinopril', 10.99, 'Dry cough, headache, fatigue', 'Take one tablet by mouth daily', '10 mg', 'ACE inhibitor', 5);
INSERT INTO medications VALUES(2, 'Atorvastatin', 15.49, 'Muscle pain, nausea, diarrhea', 'Take one tablet by mouth daily at bedtime', '20 mg', 'Statin', 4);
INSERT INTO medications VALUES(3, 'Metformin', 5.99, 'Stomach upset, diarrhea, metallic taste', 'Take one tablet by mouth twice daily with meals', '500 mg', 'Antidiabetic', 2);
INSERT INTO medications VALUES(4, 'Levothyroxine', 8.99, 'Hair loss, weight loss, sweating', 'Take one tablet by mouth daily on an empty stomach', '50 mg', 'Thyroid hormone', 3);
INSERT INTO medications VALUES(5, 'Ibuprofen', 3.49, 'Stomach upset, nausea, heartburn', 'Take one tablet by mouth every 4-6 hours as needed', '200 mg', 'Nonsteroidal anti-inflammatory drug', 1);



-- simple queries into each table
-- show the details of each doctor who has office in CA and their specialization is Dermatologist
Select * from doctors
where state = 'CA' AND specialization = 'Dermatologist'
Order by last_name;


-- listing patients older than 50 by using `months_between`
select * from patients
where trunc(months_between(sysdate, birthdate) / 12) > 50
Order by birthdate desc;


-- showing all appointments that are scheduled, and ordering by how long the appointment is
select * from appointments
where status in 'scheduled'
order by duration;


-- listing all medications where there is some side effect that relates to stomach issues
select * from medications
where LOWER(side_effects) like '%stomach%'
order by price desc;

-- Aggregate function
-- aggregate function that creates specializationCount that grabs the amount of doctors with the given specialization (Dermatologist)
select count(*) as specializationCount
from doctors
where specialization = 'Dermatologist';


-- aggregate function that creates PatientAbove50 that grabs the amount of patients older than 50
select count(*) as PatientAbove50
From patients
where Extract(Year from sysdate) - Extract(Year from birthdate) >50;

-- creates LongestAppointment, which is the appointment with the longest expected duration, and TotalDuration, which is the sum of all appointment durations
select MAX(duration) AS LongestAppointment, SUM(duration) as TotalDuration
from appointments;

-- creates MinimumPrice, which retrieves the price of the cheapest medicine, and MaxDosage, which retrieves the largest dosage (in mg)
SELECT MIN(price) as MinimumPrice, Max(dosage) as MaxDosage
from medications;


-- Insert more doctors
INSERT INTO doctors VALUES (6, 'Alice', 'Wong', 'alicewong@gmail.com', '555-666-7777', 'Orthopedist', '6789 Cedar St', NULL, 'Seattle', 'WA', '98101');
INSERT INTO doctors VALUES (7, 'Michael', 'Garcia', 'mgarcia@gmail.com', '555-888-9999', 'Gastroenterologist', '1234 Birch St', NULL, 'Miami', 'FL', '33101');


-- Insert more patients with relationships to doctors
INSERT INTO patients VALUES (6, 'Samantha', 'Adams', 'samadams@gmail.com', '555-111-2222', 'female', TO_DATE('1980-06-06', 'YYYY-MM-DD'), '9876 Willow St', NULL, 'Boston', 'MA', '02101', 6);
INSERT INTO patients VALUES (7, 'James', 'Rodriguez', 'jrodriguez@gmail.com', '555-333-4444', 'male', TO_DATE('1970-07-07', 'YYYY-MM-DD'), '4321 Spruce St', NULL, 'Denver', 'CO', '80201', 7);
INSERT INTO patients VALUES (8, 'Nancy', 'Martinez', 'nmartinez@gmail.com', '555-555-6666', 'female', TO_DATE('1960-08-08', 'YYYY-MM-DD'), '5678 Pine St', NULL, 'Atlanta', 'GA', '30301', 7);


-- Insert more appointments with relationships to doctors and patients
INSERT INTO appointments VALUES (6, 6, 6, TIMESTAMP '2023-05-06 09:00:00', 45, 'Knee pain', 'Naproxen', 'completed');
INSERT INTO appointments VALUES (7, 7, 7, TIMESTAMP '2023-05-07 10:30:00', 30, 'Stomach pain', 'Omeprazole', 'completed');
INSERT INTO appointments VALUES (8, 7, 8, TIMESTAMP '2023-05-08 11:00:00', 30, 'Heartburn', 'Ranitidine', 'scheduled');


-- Insert more medications with relationships to patients
INSERT INTO medications VALUES (6, 'Naproxen', 6.99, 'Stomach upset, nausea, headache', 'Take one tablet by mouth every 8-12 hours as needed', '500 mg', 'Nonsteroidal anti-inflammatory drug', 6);
INSERT INTO medications VALUES (7, 'Omeprazole', 12.49, 'Headache, diarrhea, stomach pain', 'Take one capsule by mouth daily before a meal', '20 mg', 'Proton pump inhibitor', 7);
INSERT INTO medications VALUES (8, 'Ranitidine', 4.99, 'Headache, dizziness, constipation', 'Take one tablet by mouth twice daily', '150 mg', 'H2 antagonist', 8);
INSERT INTO appointments VALUES (9, 2, 1, TIMESTAMP '2023-05-09 14:00:00', 30, 'Skin rash', 'Hydrocortisone', 'scheduled');
INSERT INTO appointments VALUES (10, 2, 2, TIMESTAMP '2023-05-10 15:00:00', 30, 'Eczema', 'Clobetasol', 'scheduled');
INSERT INTO appointments VALUES (11, 2, 3, TIMESTAMP '2023-05-11 16:00:00', 30, 'Psoriasis', 'Calcipotriene', 'scheduled');



-- grouping together the number of total appointments by each doctor, so we can see who currently has scheduled how many appointments
SELECT doctor_id, COUNT(*) as num_appointments
FROM appointments
GROUP BY doctor_id
ORDER BY doctor_id;

-- grabbing category, side_effects, and the amount of medications that match the category and side_effects description
SELECT category, side_effects, COUNT(*) as num_medications
FROM medications
GROUP BY category, side_effects
ORDER BY category, num_medications DESC;



-- selecting all patients that are seeing a cardiologist
SELECT *
FROM patients
WHERE id IN (
   SELECT a.patient_id
   FROM appointments a
   JOIN doctors d ON a.doctor_id = d.id
   WHERE lower(d.specialization) = 'cardiologist'
);



-- selecting from medications medicines that have a price above the average price
SELECT *
FROM medications
WHERE price > (
   SELECT AVG(price)
   FROM medications
);


-- selecting the doctor that has the earliest upcoming appointment (only the first)
SELECT *
FROM doctors
WHERE ID =
   (SELECT doctor_id
    FROM appointments
    ORDER BY date_time desc
    FETCH FIRST 1 ROWS ONLY);



-- using join to determine which doctors have appointments
SELECT first_name, last_name, specialization, patient_id, duration, date_time
FROM doctors
JOIN appointments ON doctors.id = appointments.id;


INSERT INTO doctors VALUES (8, 'John', 'Doe', 'johndoe@example.com', '555-123-4567', 'Cardiologist', '5678 Maple Ave', NULL, 'New York', 'NY', '10001');


-- using left join to help indicate when doctors do (not) have appointments
SELECT d.id, d.first_name, d.last_name, COUNT(a.id) as num_appointments
FROM doctors d
LEFT JOIN appointments a ON a.doctor_id = d.id
GROUP BY d.id, d.first_name, d.last_name;


-- Insert new values here
INSERT INTO patients VALUES(9, 'James', 'Smith', 'james.smith@example.com', '555-999-8888', 'male', TO_DATE('1980-08-08', 'YYYY-MM-DD'), '6789 Oak St', NULL, 'Los Angeles', 'CA', '90001', 5);



-- using left and right join to allow us to see if any patients currently do not have any appointments
SELECT a.id AS appointment_id, a.date_time, a.duration, d.first_name AS doctor_first_name, d.last_name AS doctor_last_name, p.first_name AS patient_first_name, p.last_name AS patient_last_name
FROM appointments a
LEFT JOIN doctors d ON a.doctor_id = d.id
RIGHT JOIN patients p ON a.patient_id = p.id;


-- PROCEDURE for cancellation
DROP PROCEDURE cancel_appointment;


CREATE OR REPLACE PROCEDURE cancel_appointment(appointment_id IN NUMBER)
IS
BEGIN
 UPDATE appointments SET status = 'canceled' WHERE id = appointment_id;
 COMMIT;
 DBMS_OUTPUT.PUT_LINE('Appointment ' || appointment_id || ' has been successfully canceled.');
END cancel_appointment;
/

-- Test the procedure
BEGIN
 cancel_appointment(1);
END;
/

-- Procedure for adding a new appointment
DROP PROCEDURE add_new_appointment;
CREATE OR REPLACE PROCEDURE add_new_appointment(
 id IN NUMBER,
 p_doctor_id IN NUMBER,
 p_patient_id IN NUMBER,
 p_date_time IN TIMESTAMP,
 p_duration IN NUMBER,
 p_reason IN VARCHAR2,
 p_prescription IN VARCHAR2,
 p_status IN VARCHAR2
)
IS
BEGIN
 INSERT INTO appointments VALUES (id, p_doctor_id, p_patient_id, p_date_time, p_duration, p_reason, p_prescription, p_status);
 COMMIT;
 DBMS_OUTPUT.PUT_LINE('Appointment ' || id || ' for doctor ' || p_doctor_id || ' and patient ' || p_patient_id || 'successfully created');
END add_new_appointment;
/

-- Test the procedure
BEGIN
 add_new_appointment(25, 1, 1, TO_TIMESTAMP('2023-05-01 15:00:00', 'YYYY-MM-DD HH24:MI:SS'), 30, 'Routine check-up', NULL, 'scheduled');
END;
/


-- Trigger to refuse addition of overlapping appointment by doctor id
CREATE OR REPLACE TRIGGER check_overlap_before_insert
BEFORE INSERT ON appointments
FOR EACH ROW
DECLARE
 overlap_count INT;
BEGIN
 SELECT COUNT(*)
 INTO overlap_count
 FROM appointments
 WHERE doctor_id = :NEW.doctor_id
   AND TRUNC(date_time, 'DD') = TRUNC(:NEW.date_time, 'DD')
   AND ((:NEW.date_time >= date_time AND :NEW.date_time < date_time + NUMTODSINTERVAL(1800, 'SECOND'))
     OR (date_time + NUMTODSINTERVAL(1800, 'SECOND') > :NEW.date_time AND date_time + NUMTODSINTERVAL(1800, 'SECOND') <= :NEW.date_time));


 IF overlap_count > 0 THEN
   raise_application_error(-20001, 'Cannot insert overlapping appointment for the same doctor.');
 END IF;
END;
/
-- Test the trigger
INSERT INTO appointments (doctor_id, patient_id, date_time, duration, reason, prescription, status)
VALUES (1, 2, TO_TIMESTAMP('2023-05-01 10:00:00', 'YYYY-MM-DD HH24:MI:SS'), 30, 'Routine check-up', NULL, 'scheduled');


-- Trigger to refuse appointment after 6pm.
CREATE OR REPLACE TRIGGER check_after_hours_before_insert
BEFORE INSERT ON appointments
FOR EACH ROW
BEGIN
 IF TO_CHAR(:NEW.date_time, 'HH24:MI:SS') >= '18:00:00' THEN
   raise_application_error(-20002, 'Cannot insert appointment after 6 PM.');
 END IF;
END;
/


-- test the trigger
INSERT INTO appointments (id, patient_id, doctor_id, date_time, status) VALUES (1, 1, 1, TO_DATE('2023-04-28 19:00:00', 'YYYY-MM-DD HH24:MI:SS'), 'scheduled');



-- Functions
-- Function to get patient id from date and time
CREATE OR REPLACE FUNCTION get_patient_id(appt_date IN DATE, appt_time IN TIMESTAMP)
RETURN NUMBER
IS
 pat_id NUMBER;
BEGIN
 SELECT patient_id INTO pat_id FROM appointments WHERE TRUNC(date_time) = TRUNC(appt_date) AND TO_CHAR(date_time, 'HH24:MI:SS') = TO_CHAR(appt_time, 'HH24:MI:SS');
 RETURN pat_id;
END;
/


-- Test the function
SELECT get_patient_id(DATE '2023-05-01', TIMESTAMP '2023-05-01 10:00:00') FROM DUAL;



-- Function to get dosage of medication from patient id and medication id
CREATE OR REPLACE FUNCTION get_dosage(
 p_pat_id IN NUMBER,
 p_med_id IN NUMBER
) RETURN VARCHAR2
IS
 dosage VARCHAR2(50);
BEGIN
 SELECT medications.dosage INTO dosage
 FROM medications
 WHERE medications.pat_id = p_pat_id AND medications.id = p_med_id;


 RETURN dosage;
END;
/


-- Test the function
SELECT get_dosage(1, 5) FROM DUAL;

-- Extra functions, triggers, and procedures.

-- Function to calculate the total duration of all appointments for a given doctor on a specific day
CREATE OR REPLACE FUNCTION get_total_duration_for_doctor(doctor_id IN NUMBER, appt_date IN DATE)
RETURN NUMBER
IS
total_duration NUMBER;
BEGIN
SELECT SUM(duration) INTO total_duration FROM appointments WHERE doctor_id = doctor_id AND TRUNC(date_time) = appt_date;
RETURN total_duration;
END;
/

-- Test the function
SELECT get_total_duration_for_doctor(1, DATE '2023-05-01') FROM DUAL;

-- Trigger to enforce a maximum number of appointments per day for a given doctor
CREATE OR REPLACE TRIGGER max_appointments_per_day
BEFORE INSERT ON appointments
FOR EACH ROW
DECLARE
appt_count NUMBER;
max_appts_per_day CONSTANT NUMBER := 5;
BEGIN
SELECT COUNT(*) INTO appt_count FROM appointments WHERE doctor_id = :NEW.doctor_id AND TRUNC(date_time) = TRUNC(:NEW.date_time);
IF appt_count >= max_appts_per_day THEN
raise_application_error(-20003, 'Maximum number of appointments per day reached for this doctor.');
END IF;
END;
/

-- Procedure to update the status of all appointments for a given patient
CREATE OR REPLACE PROCEDURE update_all_appointment_statuses_for_patient(patient_id IN NUMBER, status IN VARCHAR2)
IS
BEGIN
UPDATE appointments SET status = status WHERE patient_id = patient_id;
COMMIT;
DBMS_OUTPUT.PUT_LINE('Status updated for all appointments of patient ' || patient_id);
END update_all_appointment_statuses_for_patient;
/

-- Test the procedure
BEGIN
update_all_appointment_statuses_for_patient(1, 'canceled');
END;
/