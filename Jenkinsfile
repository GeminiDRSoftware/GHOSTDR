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

        stage ("Tests") {
            environment {
                MPLBACKEND = "agg"
                PATH = "$JENKINS_CONDA_HOME/bin:$PATH"
                DRAGONS_TEST_OUT = "./integ_tests_outputs/"
                TOX_ARGS = ""
                TMPDIR = "${env.WORKSPACE}/.tmp/integ/"
            }
            steps {
                echo "Running build #${env.BUILD_ID} on ${env.NODE_NAME}"
                checkout scm
                echo "${env.PATH}"
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