import React, { useState } from "react";
import { StandaloneStructServiceProvider } from "ketcher-standalone";
import { Editor } from "ketcher-react";
import Button from 'react-bootstrap/Button';

import classes from "./ChemEditor.module.css";
import { Container } from "react-bootstrap";

const ChemEditor = () => {
  const [smileString, setSmileString] = useState("");

  const structServiceProvider = new StandaloneStructServiceProvider();

  const handleGetSmile = () => {
    window.ketcher.getSmiles().then((smile) => {
      setSmileString(smile);
    });
  };

  return (
    <div className={classes.wrapper}>
      <div className={classes.chemeditor}>
        <Editor
          staticResourcesUrl={process.env.PUBLIC_URL}
          structServiceProvider={structServiceProvider}
          onInit={(ketcher) => {
            window.ketcher = ketcher;
            window.parent.postMessage(
              {
                eventType: "init",
              },
              "*"
            );
          }}
        />
      </div>
      <div className={classes.controls}>
        <Container>
          <Button variant="primary" onClick={handleGetSmile} className="me-2">
            Search structure
          </Button>
          <Button variant="success" onClick={handleGetSmile} className="ms-2">
            Save structure
          </Button>
        </Container>
        <div>{smileString}</div>
      </div>
    </div>
  );
};

export default ChemEditor;
