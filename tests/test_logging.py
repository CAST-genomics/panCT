import logging

from pytest import LogCaptureFixture

from panct.logging import getLogger


def test_getLogger(caplog: LogCaptureFixture):
    logger = getLogger()
    assert isinstance(logger, logging.Logger)

    # capture any messages that get logged to the Logger
    with caplog.at_level(logging.ERROR):
        logger.error("Error message")

    # check that we properly logged the message
    assert "Error message" in caplog.text


def test_getLogger_instance():
    logger = getLogger()
    assert isinstance(logger, logging.Logger)


def test_getLogger_with_name():
    logger_name = "test_logger"
    logger = getLogger(name=logger_name)
    assert logger.name == "panct.test_logger"


def test_getLogger_with_level():
    logger = getLogger(level="DEBUG")
    assert logger.level == logging.DEBUG
    logger = getLogger(level=logging.DEBUG)
    assert logger.level == logging.DEBUG


def test_getLogger_with_exact_time():
    logger = getLogger()
    formatter = logger.handlers[0].formatter
    assert "%(msecs)" not in formatter._fmt


def test_getLogger_with_index_name():
    log = getLogger(name="index")
    log.setLevel(logging.ERROR)
    assert isinstance(log, logging.Logger)
    assert log.name == "panct.index"


def test_getLogger_with_default_parameters():
    logger = getLogger()
    assert isinstance(logger, logging.Logger)
    assert logger.name == "panct"
    # by default, the log level should be "ERROR"
    assert logger.level == logging.ERROR
