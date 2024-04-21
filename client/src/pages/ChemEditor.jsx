import React, {useState} from "react";
import {StandaloneStructServiceProvider} from "ketcher-standalone";
import {Editor} from "ketcher-react";


const ChemEditor = () => {

  const [smileString, setSmileString] = useState('')

  const structServiceProvider = new StandaloneStructServiceProvider()

  const handleGetSmile = () => {
    window.ketcher.getSmiles().then(
      (smile) => {setSmileString(smile)}
    )
  }

  return (
    <div style={{height: "100%"}}>
      <div style={{height: "35%"}}>
        <Editor
          staticResourcesUrl={process.env.PUBLIC_URL}
          structServiceProvider={structServiceProvider}
          onInit={(ketcher) => {
            window.ketcher = ketcher
            window.parent.postMessage(
            {
              eventType: 'init',
            },
            '*',
          );
          }}
        />
      </div>
      <button
        onClick={handleGetSmile}
      >
        get smile
      </button>
      <div>
        {smileString}
      </div>
    </div>

  );
};

export default ChemEditor;
