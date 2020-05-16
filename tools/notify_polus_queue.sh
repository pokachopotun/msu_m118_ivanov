#!/bin/bash

jobs=$(./polus_jobs.sh)

notify-send --urgency=critical "Polus status" "JOBS = $jobs"
