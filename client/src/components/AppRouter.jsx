import React from 'react'
import { Navigate, Route, Routes } from "react-router-dom";
import ChemEditor from '../pages/chemeditor/ChemEditor';
import { Login } from '../pages/login/Login';

export const AppRouter = () => {
  return (
    <Routes>
      <Route path="/login" element={<Login />} />
      <Route path="/editor" element={<ChemEditor />} />
      <Route path="/*" element={<Navigate to="/editor" />} />
    </Routes>
  )
}
