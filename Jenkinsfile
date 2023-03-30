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

        stage ("New tests") {
            environment {
                MPLBACKEND = "agg"
                PATH = "$JENKINS_CONDA_HOME/bin:$PATH"
                DRAGONS_TEST_OUT = "./new_tests_outputs/"
                TOX_ARGS = "ghost_instruments ghostdr"
                TMPDIR = "${env.WORKSPACE}/.tmp/new/"
            }
            steps {
                echo "Running build #${env.BUILD_ID} on ${env.NODE_NAME}"
                checkout scm
                echo "${env.PATH}"
                sh '.jenkins/scripts/setup_agent.sh'
                sh 'tox -e ghost-ghostnew -v -r -- --basetemp=${DRAGONS_TEST_OUT} ${TOX_ARGS}'
            }
            post {
                always {
                    echo "Delete temporary folder: ${TMPDIR}"
                    dir ( '$TMPDIR' ) { deleteDir() }
                    echo "Delete Tox Environment: .tox/ghost-ghostnew"
                    dir ( ".tox/ghost-ghostnew" ) { deleteDir() }
                    echo "Delete Outputs folder: "
                    dir ( "${DRAGONS_TEST_OUT}" ) { deleteDir() }
                }
            }
        }

        /*
        stage ("GHOST parallel tests") {

            parallel {

                stage ("Unit tests") {
                    environment {
                        MPLBACKEND = "agg"
                        PATH = "$JENKINS_CONDA_HOME/bin:$PATH"
                        DRAGONS_TEST_OUT = "./unit_tests_outputs/"
                        TOX_ARGS = "ghost_instruments ghostdr"
                        TMPDIR = "${env.WORKSPACE}/.tmp/unit/"
                    }
                    steps {
                        echo "Running build #${env.BUILD_ID} on ${env.NODE_NAME}"
                        checkout scm
                        echo "${env.PATH}"
                        sh '.jenkins/scripts/setup_agent.sh'
                        sh 'tox -e ghost-ghostunit -v -r -- --basetemp=${DRAGONS_TEST_OUT} ${TOX_ARGS}'
                    }
                    post {
                        always {
                            echo "Delete temporary folder: ${TMPDIR}"
                            dir ( '$TMPDIR' ) { deleteDir() }
                            echo "Delete Tox Environment: .tox/ghost-ghostunit"
                            dir ( ".tox/ghost-ghostunit" ) { deleteDir() }
                            echo "Delete Outputs folder: "
                            dir ( "${DRAGONS_TEST_OUT}" ) { deleteDir() }
                        }
                    }
                }

                stage ("Bundle tests") {
                    environment {
                        MPLBACKEND = "agg"
                        PATH = "$JENKINS_CONDA_HOME/bin:$PATH"
                        DRAGONS_TEST_OUT = "./bundle_tests_outputs/"
                        TOX_ARGS = "ghost_instruments ghostdr"
                        TMPDIR = "${env.WORKSPACE}/.tmp/bundle/"
                    }
                    steps {
                        echo "Running build #${env.BUILD_ID} on ${env.NODE_NAME}"
                        sleep 5  // needed to avoid timeout from simultaneous checkout
                        checkout scm
                        echo "${env.PATH}"
                        sh '.jenkins/scripts/setup_agent.sh'
                        sh 'tox -e ghost-ghostbundle -v -r -- --basetemp=${DRAGONS_TEST_OUT} ${TOX_ARGS}'
                    }
                    post {
                        always {
                            echo "Delete temporary folder: ${TMPDIR}"
                            dir ( '$TMPDIR' ) { deleteDir() }
                            echo "Delete Tox Environment: .tox/ghost-ghostbundle"
                            dir ( ".tox/ghost-ghostbundle" ) { deleteDir() }
                            echo "Delete Outputs folder: "
                            dir ( "${DRAGONS_TEST_OUT}" ) { deleteDir() }
                        }
                    }
                }

                stage ("SLITV tests") {
                    environment {
                        MPLBACKEND = "agg"
                        PATH = "$JENKINS_CONDA_HOME/bin:$PATH"
                        DRAGONS_TEST_OUT = "./slit_tests_outputs/"
                        TOX_ARGS = "ghost_instruments ghostdr"
                        TMPDIR = "${env.WORKSPACE}/.tmp/slit/"
                    }
                    steps {
                        echo "Running build #${env.BUILD_ID} on ${env.NODE_NAME}"
                        sleep 10  // needed to avoid timeout from simultaneous checkout
                        checkout scm
                        echo "${env.PATH}"
                        sh '.jenkins/scripts/setup_agent.sh'
                        sh 'tox -e ghost-ghostslit -v -r -- --basetemp=${DRAGONS_TEST_OUT} ${TOX_ARGS}'
                    }
                    post {
                        always {
                            echo "Delete temporary folder: ${TMPDIR}"
                            dir ( '$TMPDIR' ) { deleteDir() }
                            echo "Delete Tox Environment: .tox/ghost-ghostslit"
                            dir ( ".tox/ghost-ghostslit" ) { deleteDir() }
                            echo "Delete Outputs folder: "
                            dir ( "${DRAGONS_TEST_OUT}" ) { deleteDir() }
                        }
                    }
                }

            }  // end parallel
        }  // end stage
        */


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