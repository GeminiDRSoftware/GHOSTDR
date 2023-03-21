#!/usr/bin/env groovy
/*
 * Jenkins Pipeline for GHOSTDR
 *
 * by Chris Simpson (adapted from BCQ's DRAGONS pipeline)
 *
 * Required Plug-ins:
 * - CloudBees File Leak Detector?
 * - Cobertura Plug-in?
 * - Warnings NG?
 */

pipeline {

    agent any

    options {
        skipDefaultCheckout(true)
        buildDiscarder(logRotator(numToKeepStr: '5'))
        timestamps()
        timeout(time: 4, unit: 'HOURS')
    }

    stages {

        stage ("Unit tests") {
            environment {
                MPLBACKEND = "agg"
                PATH = "$JENKINS_CONDA_HOME/bin:$PATH"
                DRAGONS_TEST_OUT = "./unit_tests_outputs/"
                TOX_ARGS = ""
                TMPDIR = "${env.WORKSPACE}/.tmp/unit/"
            }
            steps {
                echo "Running build #${env.BUILD_ID} on ${env.NODE_NAME}"
                checkout scm
                echo "${env.PATH}"
                /* sh '.jenkins/scripts/setup_agent.sh' */
                echo "Running tests with Python 3.7"
                sh 'tox -e ghost-unit -v -r -- --basetemp=${DRAGONS_TEST_OUT} --junit-xml reports/unittests_results.xml ${TOX_ARGS}'
                echo "Reportint coverage to CodeCov"
                sh 'tox -e codecov -- -F unit'
            }



        }


    }

    post {
        success {
            deleteDir() /* clean up our workspace */
        }
        failure {
            deleteDir()
        }
    }


}